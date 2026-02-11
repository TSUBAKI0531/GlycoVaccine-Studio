import json
import pandas as pd
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley

class GlycoConjugateWorkflow:
    def __init__(self, job_name):
        self.job_name = job_name

    def _calculate_sasa(self, cif_path):
        parser = MMCIFParser(QUIET=True)
        struct = parser.get_structure("model", cif_path)[0]
        sr = ShrakeRupley()
        sr.compute(struct, level="R")
        glycan_sasa = sum(getattr(res, 'sasa', 0.0) for chain in struct if chain.id != 'A' for res in chain)
        return {"glycan_sasa": glycan_sasa}

    def create_full_complex_json(self, job_name, antigen_prot, glycan_data, h_seq, l_seq, mode="Web"):
        """
        AlphaFold 3用JSON生成
        mode="Web": AlphaFold Server (Web) 用。CCDコードを使用。
        mode="Standalone": ローカル版用。SMILESを使用。
        """
        antigen_prot = antigen_prot.strip().upper()
        h_seq = h_seq.strip().upper()
        l_seq = l_seq.strip().upper()

        if mode == "Web":
            # Web Server版: CCDコードを使用し、sequences/proteinChainキーを用いる
            job_data = {
                "name": job_name,
                "modelSeeds": [1],
                "sequences": [
                    {"proteinChain": {"sequence": antigen_prot, "count": 1}},
                    {"ligand": {"ligand": glycan_data, "count": 1}}, # CCDコードを入力
                    {"proteinChain": {"sequence": h_seq, "count": 1}},
                    {"proteinChain": {"sequence": l_seq, "count": 1}}
                ]
            }
        else:
            # Standalone版: SMILESを使用し、modelContentsキーを用いる
            job_data = {
                "name": job_name,
                "modelContents": [
                    {"protein": {"sequence": antigen_prot}},
                    {"ligand": {"smiles": glycan_data}}, # SMILESを入力
                    {"protein": {"sequence": h_seq}},
                    {"protein": {"sequence": l_seq}}
                ]
            }
        return [job_data] if mode == "Web" else job_data

# --- 他のクラス(CDRScorer, AntibodyDesigner等)は維持 ---