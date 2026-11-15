import csv
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import sys
import types

from rnaends2tracks.conf import project_from_conf, ConfError
from rnaends2tracks.config import workflow_requirements
from rnaends2tracks.gex_annotation import build_annotation
from rnaends2tracks.gex import export_counts
from rnaends2tracks.enrichment import _extended_gene_table
from rnaends2tracks.gex_report import enrichment_summary
from rnaends2tracks.alignment_inputs import validate_reused_alignments


def read(path):
    with path.open(newline="",encoding="utf-8") as handle:
        return list(csv.DictReader(handle,delimiter="\t"))


class ExtendedExpressionTests(unittest.TestCase):
    def test_annotation_both_strands_promoters_nested_genes_and_original_exons(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gtf = root / "genes.gtf"
            records = []
            def gene(name, start, end, strand, exons=None, biotype="protein_coding", tss=None):
                for feature,left,right in [("gene",start,end),("transcript",tss or start,end),
                                           *[("exon",a,b) for a,b in (exons or [(start,end)])]]:
                    records.append(f'chr1\ttest\t{feature}\t{left}\t{right}\t.\t{strand}\t.\tgene_id "{name}"; gene_type "{biotype}";\n')
            gene("plus",101,200,"+",[(101,150),(140,200)])
            gene("next_plus",3201,3300,"+")
            gene("opposite",500,550,"-")
            gene("minus",8001,8100,"-")
            gene("next_minus",4901,5000,"-")
            gene("nested",8150,8200,"+", tss=8150)
            gene("outer",7900,9000,"+",[(7900,7950),(8950,9000)])
            gene("edge",9901,10000,"+")
            gene("smallRNA",9501,9520,"+",biotype="miRNA")
            gtf.write_text("".join(records))
            sizes = root / "sizes"; sizes.write_text("chr1\t10000\n")
            outputs = build_annotation(gtf,sizes,root / "annotation")
            audit = {r["gene_id"]:r for r in read(outputs["audit"])}
            self.assertEqual(audit["plus"]["extension_end"],"1201")
            self.assertEqual(audit["plus"]["limiting_gene"],"next_plus")
            self.assertIn("opposite",audit["plus"]["opposite_strand_extension_overlaps"])
            self.assertEqual(audit["minus"]["extension_start"],"7000")
            self.assertEqual(audit["nested"]["accepted_extension_nt"],"0")
            self.assertEqual(audit["edge"]["accepted_extension_nt"],"0")
            self.assertEqual(audit["smallRNA"]["accepted_extension_nt"],"0")
            saf = read(outputs["saf"])
            plus = [r for r in saf if r["GeneID"] == "plus"]
            self.assertEqual([(r["Start"],r["End"]) for r in plus],[("101","1201")])
            before = outputs["saf"].read_bytes()
            build_annotation(gtf,sizes,root / "annotation")
            self.assertEqual(outputs["saf"].read_bytes(),before)
            gtf.write_text("".join(reversed(records)))
            reordered = build_annotation(gtf,sizes,root / "reordered")
            self.assertEqual(reordered["saf"].read_bytes(),before)
            outputs["saf"].write_text("corrupt")
            build_annotation(gtf,sizes,root / "annotation")
            self.assertEqual(outputs["saf"].read_bytes(),before)

    def test_default_switches_and_bam_input_requirements(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = Path(tmp) / "config.conf"
            base = "PROJECT_ID=x\nSAMPLESHEET=samples.csv\nOUTPUT_DIR=results\nRUN_APA_A_MCELL2019=false\nRUN_APA_A2=false\nRUN_TRACKS=false\n"
            config.write_text(base)
            project,_ = project_from_conf(config)
            self.assertEqual(project["gene_expression"]["methods"],["extended_exons_v1"])
            self.assertTrue(project["gene_expression"]["stringent"])
            self.assertTrue(project["gene_expression"]["relaxed"])
            self.assertFalse(workflow_requirements(project)["active_pas"])
            config.write_text(base+"INPUT_MODE=bam\nALIGNMENT_SOURCE_DIR=old/02_alignment\nCLEANUP_INTERMEDIATES=false\nRUN_RSEQC=false\n")
            project,_ = project_from_conf(config)
            self.assertEqual(project["input_mode"],"bam")
            config.write_text(base+"DGE_METHODS=extended_exons_v1,legacy_c4\n")
            project,_ = project_from_conf(config)
            self.assertTrue(workflow_requirements(project)["active_pas"])
            config.write_text(base+"DGE_METHODS=legacy_c4\n")
            with self.assertRaises(ConfError):
                project_from_conf(config)  # primary must also be selected

    def test_featurecounts_sample_identity_and_integral_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); bams = [root / "A.bam",root / "B.bam"]
            source,target = root / "featureCounts.tsv",root / "counts.tsv"
            header = "Geneid\tChr\tStart\tEnd\tStrand\tLength\t"+"\t".join(map(str,bams))+"\n"
            source.write_text(header+"g1\tchr1\t1\t10\t+\t10\t4\t8\n")
            export_counts(source,target,["A","B"],bams)
            self.assertEqual(read(target),[{"gene_id":"g1","A":"4","B":"8"}])
            with self.assertRaises(RuntimeError):
                export_counts(source,target,["B","A"],bams[::-1])
            source.write_text(header+"g1\tchr1\t1\t10\t+\t10\t4.5\t8\n")
            with self.assertRaises(ValueError):
                export_counts(source,target,["A","B"],bams)

    def test_enrichment_modes_use_formal_padj_and_common_valid_universe(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); source=root / "result.tsv"; target=root / "prepared.tsv"
            source.write_text("gene_id\tlog2FoldChange\tstat\tpadj\tvalid_test\nsmall\t0.4\t3\t0.01\tTRUE\nbig\t2\t1\t0.2\tTRUE\nbad\t3\tNA\tNA\tFALSE\n")
            self.assertEqual(_extended_gene_table(source,target,.05,"relaxed_fdr"),(2,1))
            self.assertEqual(read(target)[0]["foreground"],"1")
            self.assertEqual(_extended_gene_table(source,target,.05,"shared_gsea"),(2,0))
            self.assertEqual([r["rank_score"] for r in read(target)],["3.0","1.0"])
            gsea=root / "gsea.tsv"
            gsea.write_text("database\tterm_id\tterm_name\tNES\tpadj\nH\tp1\tPositive1\t3\t0.001\nH\tp2\tPositive2\t2\t0.002\nH\tn1\tNegative\t-2\t0.04\n")
            summary=enrichment_summary({"gsea_file":str(gsea),"numerator":"A","denominator":"B"},"GSEA",.05,1)
            self.assertEqual({r["term_id"] for r in summary},{"p1","n1"})
            self.assertIn("B",next(r for r in summary if r["term_id"]=="n1")["direction"])

    def test_reuse_input_cannot_be_output_or_nested(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                validate_reused_alignments({"alignment_source_dir":tmp,"output_dir":str(Path(tmp)/"new")},[],{},False)

    def test_bam_validation_missing_header_and_c0_policy_guards(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); source=root/"old"; output=root/"new"
            project={"input_mode":"bam","alignment_source_dir":str(source),"output_dir":str(output)}
            sample={"sample_id":"s1","genome":"GRCm39"}
            sizes=root/"sizes"; sizes.write_text("chr1\t1000\n")
            refs={"GRCm39":{"chrom_sizes":str(sizes)}}
            read_record=types.SimpleNamespace(is_unmapped=False,is_secondary=False,is_supplementary=False,
                has_tag=lambda key: True,get_tag=lambda key: 1)
            header={"HD":{"SO":"coordinate"},"RG":[{"SM":"s1"}]}
            class Bam:
                references=["chr1"]; lengths=[1000]
                def __init__(self,*args,**kwargs): self.header=types.SimpleNamespace(to_dict=lambda: header)
                def __enter__(self): return self
                def __exit__(self,*args): pass
                def fetch(self,**kwargs): return iter([read_record])
            stub=types.SimpleNamespace(quickcheck=lambda path:None,AlignmentFile=Bam)
            with patch.dict(sys.modules,{"pysam":stub}):
                with self.assertRaisesRegex(ValueError,"Missing source BAM"):
                    validate_reused_alignments(project,[sample],refs)
                bam=source/"s1"/"s1.bam"; bam.parent.mkdir(parents=True); bam.write_bytes(b"fixture")
                Path(str(bam)+".bai").write_bytes(b"index")
                self.assertEqual(validate_reused_alignments(project,[sample],refs)[0]["sampled_c0_records"],1)
                read_record.is_supplementary=True
                with self.assertRaisesRegex(ValueError,"C0 policy"):
                    validate_reused_alignments(project,[sample],refs)
                read_record.is_supplementary=False
                header["RG"][0]["SM"]="other"
                with self.assertRaisesRegex(ValueError,"sample identity"):
                    validate_reused_alignments(project,[sample],refs)


if __name__ == "__main__":
    unittest.main()
