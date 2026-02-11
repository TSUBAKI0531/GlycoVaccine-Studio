import numpy as np

class ComplexBuilder:
    """CueMol2 で Ribbon 表示が可能な PDB を生成する"""
    
    def _generate_proper_coords(self, sequence, chain_id="A"):
        """
        CueMol2 の SplineRenderer が認識可能な連続した座標を生成する。
        N -> CA -> C の順に Z 軸方向に進むように配置します。
        """
        lines = []
        atom_idx = 1
        
        # 螺旋のパラメータ (曲率を与えて Ribbon 演算を可能にする)
        radius = 2.0
        pitch = 3.0   # 残基ごとの上昇量 (CA-CA距離に近似)
        theta_step = 60 * (np.pi / 180) # 60度ずつ回転

        for i, aa in enumerate(sequence):
            res_idx = i + 1
            theta = i * theta_step
            
            # 基準となる CA 座標
            z_ca = i * pitch
            x_ca = radius * np.cos(theta)
            y_ca = radius * np.sin(theta)
            
            # 主鎖原子の相対配置 (N, CA, C, O) 
            # 鎖が Z 軸方向に一方向に進むように Z オフセットを調整
            atoms = [
                ("N ", x_ca - 0.4, y_ca + 1.0, z_ca - 1.2, "N"),
                ("CA", x_ca,       y_ca,       z_ca,       "C"),
                ("C ", x_ca + 0.4, y_ca - 1.0, z_ca + 1.2, "C"),
                ("O ", x_ca + 1.4, y_ca - 1.4, z_ca + 1.2, "O")
            ]
            
            for name, x, y, z, elem in atoms:
                # PDB ATOM レコードの厳格なカラムフォーマット
                line = f"ATOM  {atom_idx:>5}  {name}  ALA {chain_id}{res_idx:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           {elem}"
                lines.append(line)
                atom_idx += 1
        return "\n".join(lines) + "\nTER"

    def build_complex_pdb(self, protein_seq, linker_smi, glycan_smi):
        """統合 PDB 生成"""
        protein_seq = protein_seq.strip().upper()
        header = [
            "REMARK   Built by GlycoVaccine Studio v3.2",
            f"REMARK   Carrier: {protein_seq[:20]}...",
            f"REMARK   Linker: {linker_smi}",
            f"REMARK   Glycan: {glycan_smi}"
        ]
        pdb_body = self._generate_proper_coords(protein_seq)
        return "\n".join(header) + "\n" + pdb_body + "\nEND"