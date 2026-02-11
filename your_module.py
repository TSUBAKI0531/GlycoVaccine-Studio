import numpy as np

class ComplexBuilder:
    """CueMol2 で Ribbon 表示が可能な PDB を生成する"""
    
    def _generate_helical_coords(self, sequence, chain_id="A", offset_z=0.0):
        """
        リボン演算が可能なように、螺旋状の座標を生成する。
        完全な直線だと CueMol2 の SplineRenderer がエラーになるため。
        """
        lines = []
        atom_idx = 1
        
        # 螺旋のパラメータ (Alpha-helix の近似値)
        radius = 2.3  # 螺旋の半径
        pitch = 1.5   # 残基ごとの Z 軸上昇量
        theta_step = 100 * (np.pi / 180) # 残基ごとの回転角 (約100度)

        for i, aa in enumerate(sequence):
            res_idx = i + 1
            theta = i * theta_step
            
            # CA 原子の座標計算
            ca_x = radius * np.cos(theta)
            ca_y = radius * np.sin(theta)
            ca_z = offset_z + (i * pitch)
            
            # 主鎖原子セット (N, CA, C, O)
            # CA を基準に少しずつずらして配置し、ペプチド平面の向きを定義する
            atoms = [
                ("N ", ca_x - 0.5, ca_y + 1.2, ca_z - 1.0, "N"),
                ("CA", ca_x,       ca_y,       ca_z,       "C"),
                ("C ", ca_x + 0.5, ca_y - 1.2, ca_z + 1.0, "C"),
                ("O ", ca_x + 1.5, ca_y - 1.5, ca_z + 1.0, "O")
            ]
            
            for name, x, y, z, elem in atoms:
                # PDB ATOM レコードの厳密なフォーマット
                line = f"ATOM  {atom_idx:>5}  {name}  ALA {chain_id}{res_idx:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           {elem}"
                lines.append(line)
                atom_idx += 1
        return "\n".join(lines) + "\nTER"

    def build_complex_pdb(self, protein_seq, linker_smi, glycan_smi):
        """統合された PDB テキストを生成"""
        protein_seq = protein_seq.strip().upper()
        header = [
            "REMARK   Built by GlycoVaccine Studio v3.1",
            f"REMARK   Carrier: {protein_seq[:20]}...",
            f"REMARK   Linker: {linker_smi}",
            f"REMARK   Glycan: {glycan_smi}"
        ]
        pdb_body = self._generate_helical_coords(protein_seq)
        return "\n".join(header) + "\n" + pdb_body + "\nEND"