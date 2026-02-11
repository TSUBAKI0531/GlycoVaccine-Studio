import numpy as np
import pandas as pd

class ComplexBuilder:
    """Tab 1: アプリ内 3D 表示用のヘリックス構造生成"""
    def _generate_coords(self, sequence, offset_z=0.0):
        lines = []
        atom_idx = 1
        radius, pitch = 2.3, 1.5
        theta_step = 100 * (np.pi / 180)

        for i, aa in enumerate(sequence):
            theta = i * theta_step
            ca_x, ca_y, ca_z = radius * np.cos(theta), radius * np.sin(theta), offset_z + (i * pitch)
            atoms = [
                (" N  ", ca_x - 0.5, ca_y + 1.2, ca_z - 1.0, "N"),
                (" CA ", ca_x,       ca_y,       ca_z,       "C"),
                (" C  ", ca_x + 0.5, ca_y - 1.2, ca_z + 1.0, "C"),
                (" O  ", ca_x + 1.5, ca_y - 1.5, ca_z + 1.0, "O")
            ]
            for name, x, y, z, elem in atoms:
                lines.append(f"ATOM  {atom_idx:>5} {name} ALA A{i+1:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           {elem}")
                atom_idx += 1
        return "\n".join(lines) + "\nTER"

    def build_pdb(self, seq, l_smi, g_smi):
        header = f"REMARK   Antigen-Glycan Complex\nREMARK   SMILES: {l_smi}.{g_smi}"
        return header + "\n" + self._generate_coords(seq) + "\nEND"

class AntibodyGraftingEngine:
    """Tab 2: トラスツズマブのフレームワークへの CDR 移植"""
    def __init__(self):
        # トラスツズマブ (Herceptin) の FR 配列
        self.H_FR = {"F1":"EVQLVESGGGLVQPGGSLRLSCAAS", "F2":"WVRQAPGKGLEWVA", "F3":"RFTISADTSKNTAYLQMNSLRAEDTAVYYC", "F4":"WGQGTLVTVSS"}
        self.L_FR = {"F1":"DIQMTQSPSSLSASVGDRVTITC", "F2":"WYQQKPGKAPKLLIY", "F3":"GVPSRFSGSGSGTDFTLTISSLQPEDFATYYC", "F4":"FGQGTKVEIK"}

    def predict_cdrs(self, glycan_smi):
        """標的糖鎖に対する CDR 配列の予測ロジック (デモ用)"""
        # Tn抗原 等に基づいた配列を返す
        h_cdrs = ["GFTFSRYT", "ISSSGGST", "ARTVRYGMDV"]
        l_cdrs = ["QSVSSY", "DAS", "QQRSSWPFT"]
        return h_cdrs, l_cdrs

    def graft(self, h_cdrs, l_cdrs):
        """予測 CDR を移植して完全な重鎖・軽鎖を生成"""
        h = f"{self.H_FR['F1']}{h_cdrs[0]}{self.H_FR['F2']}{h_cdrs[1]}{self.H_FR['F3']}{h_cdrs[2]}{self.H_FR['F4']}"
        l = f"{self.L_FR['F1']}{l_cdrs[0]}{self.L_FR['F2']}{l_cdrs[1]}{self.L_FR['F3']}{l_cdrs[2]}{self.L_FR['F4']}"
        return h, l