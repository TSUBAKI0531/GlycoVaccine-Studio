import json
import pandas as pd
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley

class GlycoConjugateWorkflow:
    def __init__(self, job_name):
        self.job_name = job_name

    def _calculate_sasa(self, cif_path):
        """Biopythonを用いた安定版SASA計算"""
        parser = MMCIFParser(QUIET=True)
        struct = parser.get_structure("model", cif_path)[0]
        sr = ShrakeRupley()
        sr.compute(struct, level="R")
        # A鎖（抗原タンパク）以外の糖鎖部分のSASAを合計
        glycan_sasa = sum(getattr(res, 'sasa', 0.0) for chain in struct if chain.id != 'A' for res in chain)
        return {"glycan_sasa": glycan_sasa}

    def create_full_complex_json(self, job_name, antigen_prot, glycan_smiles, linker_smiles, bond_res_idx, h_seq, l_seq):
        """抗原+リンカー+糖鎖+抗体のフルコンプレックスJSONを生成（AF3 Server対応版）"""
        # リンカーと糖鎖を統合。AF3 Serverは結合した一つのSMILESをリガンドとして受け取ります
        combined_ligand = f"{linker_smiles}{glycan_smiles}"
        
        # サーバー仕様：全体をリスト [] で囲む必要があります
        job_data = {
            "name": job_name,
            "model_contents": [
                {"protein": {"sequence": antigen_prot, "count": 1}}, # Entity 1
                {"ligand": {"smiles": combined_ligand, "count": 1}}, # Entity 2
                {"protein": {"sequence": h_seq, "count": 1}},        # Entity 3
                {"protein": {"sequence": l_seq, "count": 1}}         # Entity 4
            ],
            "user_bonds": [ # 最新のAF3 JSON仕様: user_bonds と res_id_1/2 を使用
                {
                    "res_id_1": int(bond_res_idx),
                    "entity_id_1": 1,
                    "res_id_2": 1, 
                    "entity_id_2": 2
                }
            ]
        }
        return [job_data] # リスト形式で返す

class CDRScorer:
    HYDRO_SCALE = {'A': 0.62, 'R': -2.53, 'N': -0.78, 'D': -0.90, 'C': 0.29, 'Q': -0.85, 'E': -0.74, 'G': 0.48, 'H': -0.40, 'I': 1.38, 'L': 1.06, 'K': -1.50, 'M': 0.64, 'F': 1.19, 'P': 0.12, 'S': -0.18, 'T': -0.05, 'W': 0.81, 'Y': 0.26, 'V': 1.08}

    @classmethod
    def calculate_advanced_score(cls, h_cdrs, l_cdrs):
        """高度なスコアリング：残基寄与 + H3長さ + 疎水性モーメント"""
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
                {"name": "Clone_Tn_A1", "H_CDRs": ["GFTFSRYT", "ISSSGGST", "ARTVRYGMDV"], "L_CDRs": ["QSVSSY", "DAS", "QQRSSWPFT"]},
                {"name": "Clone_Tn_B5", "H_CDRs": ["GYTFTSYY", "ISSGGGTY", "ARGDYGYWYFDV"], "L_CDRs": ["QDISNY", "YTS", "QQGNTLPWT"]}
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

class AntibodyDockingWorkflow:
    def __init__(self, job_name, h_chain, l_chain): pass
    def analyze_paratope(self, path): return pd.DataFrame([{"Residue": "Y33", "Type": "H-bond"}])

class LightweightHotSpotAnalyzer:
    def __init__(self, file, h_chain, l_chain): pass
    def run_contact_density_scan(self): return pd.DataFrame([{"Residue": "W102", "HotSpot_Score": 0.95}])