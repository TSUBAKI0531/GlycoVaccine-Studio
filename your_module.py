import json
import pandas as pd
import numpy as np

class ComplexBuilder:
    """Tab 1, 3: CueMol2 がリボン表示できる PDB 形式の作製"""
    
    def _generate_backbone_atoms(self, sequence, chain_id, offset_z=0.0):
        """CueMol2 の SplineRenderer 用に N, CA, C 原子を生成"""
        lines = []
        atom_count = 1
        for i, aa in enumerate(sequence):
            res_num = i + 1
            z_base = offset_z + (i * 3.8) # 残基間隔 3.8A
            
            # 各残基に N, CA, C の 3 原子を配置 (簡易的な直線構造)
            atoms = [
                ("N ", 0.0, 1.4, z_base - 1.0, "N"),
                ("CA", 0.0, 0.0, z_base,       "C"),
                ("C ", 0.0, 1.4, z_base + 1.0, "C")
            ]
            for atom_name, x, y, z, element in atoms:
                line = f"ATOM  {atom_count:>5}  {atom_name}  ALA {chain_id}{res_num:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           {element}"
                lines.append(line)
                atom_count += 1
        return lines

    def build_antigen_pdb(self, prot_seq, linker_smi, glycan_smi):
        """抗原複合体の PDB 生成"""
        header = f"REMARK   Built by GlycoVaccine Studio\nREMARK   SMILES: {linker_smi}.{glycan_smi}"
        coords = self._generate_backbone_atoms(prot_seq, "A")
        return header + "\n" + "\n".join(coords) + "\nTER\nEND"

    def merge_for_cuemol(self, ant_seq, h_seq, l_seq):
        """CueMol2 用の統合 PDB (抗原+抗体重鎖+軽鎖)"""
        lines = ["REMARK   Merged Complex for CueMol2 Visualization"]
        lines.extend(self._generate_backbone_atoms(ant_seq, "A"))
        lines.append("TER")
        lines.extend(self._generate_backbone_atoms(h_seq, "H", offset_z=100.0))
        lines.append("TER")
        lines.extend(self._generate_backbone_atoms(l_seq, "L", offset_z=200.0))
        lines.append("TER")
        lines.append("END")
        return "\n".join(lines)

class GlycoConjugateWorkflow:
    def create_full_complex_json(self, job_name, antigen_prot, glycan_smiles, linker_smiles, bond_res_idx, h_seq, l_seq, mode="Web"):
        """AlphaFold Server 完全準拠版 JSON (空白除去ロジック搭載)"""
        # 末尾の空白を確実に削除
        antigen_prot = antigen_prot.strip().upper()
        h_seq = h_seq.strip().upper()
        l_seq = l_seq.strip().upper()
        combined_ligand = (linker_smiles.strip() + "." + glycan_smiles.strip()).strip(".")

        if "Web" in mode:
            job_data = {
                "name": job_name,
                "modelSeeds": [1],
                "sequences": [
                    {"proteinChain": {"sequence": antigen_prot, "count": 1}},
                    {"ligand": {"smiles": combined_ligand, "count": 1}},
                    {"proteinChain": {"sequence": h_seq, "count": 1}},
                    {"proteinChain": {"sequence": l_seq, "count": 1}}
                ]
            }
            return [job_data]
        else:
            job_data = {
                "name": job_name,
                "modelContents": [
                    {"protein": {"sequence": antigen_prot, "label": "Antigen"}},
                    {"ligand": {"smiles": combined_ligand, "label": "Glycan"}},
                    {"protein": {"sequence": h_seq, "label": "H-chain"}},
                    {"protein": {"sequence": l_seq, "label": "L-chain"}}
                ]
            }
            return job_data

class AntibodyGraftingEngine:
    def __init__(self):
        self.H_FR = {"F1":"EVQLVESGGGLVQPGGSLRLSCAAS", "F2":"WVRQAPGKGLEWVA", "F3":"RFTISADTSKNTAYLQMNSLRAEDTAVYYC", "F4":"WGQGTLVTVSS"}
        self.L_FR = {"F1":"DIQMTQSPSSLSASVGDRVTITC", "F2":"WYQQKPGKAPKLLIY", "F3":"GVPSRFSGSGSGTDFTLTISSLQPEDFATYYC", "F4":"FGQGTKVEIK"}
    def graft(self, h_cdrs, l_cdrs):
        h = f"{self.H_FR['F1']}{h_cdrs[0]}{self.H_FR['F2']}{h_cdrs[1]}{self.H_FR['F3']}{h_cdrs[2]}{self.H_FR['F4']}"
        l = f"{self.L_FR['F1']}{l_cdrs[0]}{self.L_FR['F2']}{l_cdrs[1]}{self.L_FR['F3']}{l_cdrs[2]}{self.L_FR['F4']}"
        return h, l

class CDRPredictor:
    def predict(self, smiles): return ["GFTFSRYT", "ISSSGGST", "ARTVRYGMDV"], ["QSVSSY", "DAS", "QQRSSWPFT"]