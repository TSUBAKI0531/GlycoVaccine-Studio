import pandas as pd

class ComplexBuilder:
    """タンパク質と糖鎖の情報を統合し、CueMol2 準拠の PDB を作成する"""
    
    def _generate_backbone_pdb(self, sequence, chain_id="A", offset_z=0.0):
        """CueMol2 のリボン表示に必要な N, CA, C, O 原子を生成する"""
        lines = []
        atom_idx = 1
        # アミノ酸配列をループして原子座標を計算
        for i, aa in enumerate(sequence):
            res_idx = i + 1
            z_base = offset_z + (i * 3.8) # 残基間隔 3.8A
            
            # 主鎖原子 (N, CA, C, O) の簡易的な幾何配置
            # これにより CueMol2 の SplineRenderer が正常に計算可能になる
            atoms = [
                ("N ", 0.0, 1.4, z_base - 1.0, "N"),
                ("CA", 0.0, 0.0, z_base,       "C"),
                ("C ", 0.0, 1.4, z_base + 1.0, "C"),
                ("O ", 1.2, 2.0, z_base + 1.0, "O")
            ]
            
            for name, x, y, z, elem in atoms:
                # PDB ATOM レコードの厳格なフォーマット (カラム位置を固定)
                line = f"ATOM  {atom_idx:>5}  {name}  ALA {chain_id}{res_idx:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           {elem}"
                lines.append(line)
                atom_idx += 1
        return "\n".join(lines) + "\nTER"

    def build_complex_pdb(self, protein_seq, linker_smi, glycan_smi):
        """抗原、リンカー、糖鎖情報を統合した PDB 文字列を生成"""
        # 不要な空白などをクリーニング
        protein_seq = protein_seq.strip().upper()
        linker_smi = linker_smi.strip()
        glycan_smi = glycan_smi.strip()
        
        # PDB ヘッダー情報の作成
        header = [
            "REMARK   Built by GlycoVaccine Studio",
            f"REMARK   Carrier Protein Length: {len(protein_seq)}",
            f"REMARK   Linker SMILES: {linker_smi}",
            f"REMARK   Glycan SMILES: {glycan_smi}"
        ]
        
        pdb_body = self._generate_backbone_pdb(protein_seq)
        return "\n".join(header) + "\n" + pdb_body + "\nEND"