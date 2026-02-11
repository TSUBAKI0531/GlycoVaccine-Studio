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

    def create_full_complex_json(self, job_name, antigen_prot, glycan_smiles, linker_smiles, bond_res_idx, h_seq, l_seq, platform="Web"):
        """
        AlphaFold 3用JSON生成（TypeError回避版）
        platform="Web": AlphaFold Server (Web) 用。sequences/proteinChainキーを使用。
        platform="Standalone": ローカル/スタンドアロン版用。modelContentsキーを使用。
        """
        antigen_prot = antigen_prot.strip().upper()
        h_seq = h_seq.strip().upper()
        l_seq = l_seq.strip().upper()
        
        # SMILESの結合処理（Web版でエラーの原因となる"."を除去または適切に処理）
        g_smi = glycan_smiles.strip()
        l_smi = linker_smiles.strip()
        combined_smiles = f"{l_smi}{g_smi}" if l_smi else g_smi

        if "Web" in platform:
            # Web Server版のスキーマ
            job_data = {
                "name": job_name,
                "modelSeeds": [],
                "sequences": [
                    {"proteinChain": {"sequence": antigen_prot, "count": 1}},
                    {"ligand": {"smiles": combined_smiles, "count": 1}},
                    {"proteinChain": {"sequence": h_seq, "count": 1}},
                    {"proteinChain": {"sequence": l_seq, "count": 1}}
                ]
            }
            return [job_data] # Web版はリスト形式
        else:
            # Standalone版のスキーマ
            job_data = {
                "name": job_name,
                "modelContents": [
                    {"protein": {"sequence": antigen_prot, "label": "Antigen"}},
                    {"ligand": {"smiles": combined_smiles, "label": "Glycan"}},
                    {"protein": {"sequence": h_seq, "label": "H-chain"}},
                    {"protein": {"sequence": l_seq, "label": "L-chain"}}
                ],
                "userBonds": [{"resId1": int(bond_res_idx), "entityId1": 1, "resId2": 1, "entityId2": 2}]
            }
            return job_data

# --- CDRScorer, AntibodyDesigner 等のクラスはそのまま維持してください ---
class CDRScorer:
    HYDRO_SCALE = {'A': 0.62, 'R': -2.53, 'N': -0.78, 'D': -0.90, 'C': 0.29, 'Q': -0.85, 'E': -0.74, 'G': 0.48, 'H': -0.40, 'I': 1.38, 'L': 1.06, 'K': -1.50, 'M': 0.64, 'F': 1.19, 'P': 0.12, 'S': -0.18, 'T': -0.05, 'W': 0.81, 'Y': 0.26, 'V': 1.08}
    @classmethod
    def calculate_advanced_score(cls, h_cdrs, l_cdrs):
        h3_seq = h_cdrs[2]
        all_cdrs = "".join(h_cdrs + l_cdrs)
        s_res = (all_cdrs.count('Y') + all_cdrs.count('W')) * 3.0 + (all_cdrs.count('S') + all_cdrs.count('T')) * 1.5
        h3_len = len(h3_seq)
        s_len = 5.0 if 10 <= h3_len <= 16 else -3.0
        h_values = [cls.HYDRO_SCALE.get(aa, 0) for aa in h3_seq]
        s_moment = 10.0 * (1.0 - abs(np.mean(h_values) - 0.1)) if h_values else 0
        return round(s_res + s_len + s_moment, 2)

class AntibodyDesigner:
    def __init__(self):
        self.FR_H = {"FR1": "EVQLVESGGGLVQPGGSLRLSCAAS", "FR2": "WVRQAPGKGLEWVA", "FR3": "RFTISADTSKNTAYLQMNSLRAEDTAVYYC", "FR4": "WGQGTLVTVSS"}
        self.FR_L = {"FR1": "DIQMTQSPSSLSASVGDRVTITC", "FR2": "WYQQKPGKAPKLLIY", "FR3": "GVPSRFSGSGSGTDFTLTISSLQPEDFATYYC", "FR4": "FGQGTKVEIK"}
        self.library = {
            "Tn Antigen": [
                {"name": "Clone_Tn_B5", "H_CDRs": ["GYTFTSYY", "ISSGGGTY", "ARGDYGYWYFDV"], "L_CDRs": ["QDISNY", "YTS", "QQGNTLPWT"]},
                {"name": "Clone_Tn_A1", "H_CDRs": ["GFTFSRYT", "ISSSGGST", "ARTVRYGMDV"], "L_CDRs": ["QSVSSY", "DAS", "QQRSSWPFT"]}
            ]
        }
    def get_ranked_candidates(self, motif):
        candidates = self.library.get(motif, [])
        scored_list = []
        for cand in candidates:
            score = CDRScorer.calculate_advanced_score(cand["H_CDRs"], cand["L_CDRs"])
            h_seq = f"{self.FR_H['FR1']}{cand['H_CDRs'][0]}{self.FR_H['FR2']}{cand['H_CDRs'][1]}{self.FR_H['FR3']}{cand['H_CDRs'][2]}{self.FR_H['FR4']}"
            l_seq = f"{self.FR_L['FR1']}{cand['L_CDRs'][0]}{self.FR_L['FR2']}{cand['L_CDRs'][1]}{self.FR_L['FR3']}{cand['L_CDRs'][2]}{self.FR_L['FR4']}"
            scored_list.append({**cand, "Score": score, "H_AA": h_seq, "L_AA": l_seq})
        return sorted(scored_list, key=lambda x: x["Score"], reverse=True)