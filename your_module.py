import pandas as pd

class ComplexBuilder:
    """タンパク質と糖鎖の情報を統合し、CueMol2 準拠の PDB を作成する"""
    
    def _generate_backbone_pdb(self, sequence, chain_id="A", offset_z=0.0):
        """リボン表示に必要な N, CA, C, O 原子を生成する"""
        lines = []
        atom_idx = 1
        for i, aa in enumerate(sequence):
            res_idx = i + 1
            z_base = offset_z + (i * 3.8) # 残基間隔を 3.8A に設定
            
            # 各残基の主鎖原子（N, CA, C, O）の相対座標
            # CueMol2 のレンダリングを安定させるための簡易的な幾何配置
            atoms = [
                ("N ", 0.0, 1.4, z_base - 1.0, "N"),
                ("CA", 0.0, 0.0, z_base,       "C"),
                ("C ", 0.0, 1.4, z_base + 1.0, "C"),
                ("O ", 1.2, 2.0, z_base + 1.0, "O")
            ]
            
            for name, x, y, z, elem in atoms:
                # PDB ATOM レコードの厳格なフォーマットに従う
                line = f"ATOM  {atom_idx:>5}  {name}  ALA {chain_id}{res_idx:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           {elem}"
                lines.append(line)
                atom_idx += 1
        return "\n".join(lines) + "\nTER"

    def build_complex_pdb(self, protein_seq, linker_smi, glycan_smi):
        """ヘッダー情報と主鎖座標を統合した PDB 文字列を生成"""
        # 末尾の空白などの不要な文字を除去
        protein_seq = protein_seq.strip().upper()
        linker_smi = linker_smi.strip()
        glycan_smi = glycan_smi.strip()
        
        header = [
            "REMARK   Built by GlycoVaccine Studio",
            f"REMARK   Carrier Protein Length: {len(protein_seq)}",
            f"REMARK   Linker SMILES: {linker_smi}",
            f"REMARK   Glycan SMILES: {glycan_smi}"
        ]
        
        pdb_body = self._generate_backbone_pdb(protein_seq)
        return "\n".join(header) + "\n" + pdb_body + "\nEND"