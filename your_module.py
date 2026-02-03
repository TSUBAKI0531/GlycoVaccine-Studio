import json
import os
import pandas as pd
import freesasa
from rdkit import Chem
from Bio.PDB import MMCIFParser, PDBParser, NeighborSearch
from Bio.SeqUtils import seq1

# ANARCIのインポート（環境にない場合はCDRラベルをスキップ）
try:
    from anarci import anarci
except ImportError:
    anarci = None

class GlycoConjugateWorkflow:
    """糖鎖抱合抗原の設計とモデル評価を担当するクラス"""
    def __init__(self, job_name):
        self.job_name = job_name

    def _get_af3_atom_name(self, mol, target_idx):
        """RDKitの原子IndexをAF3形式（元素記号+出現順）に変換"""
        atom = mol.GetAtomWithIdx(target_idx)
        symbol = atom.GetSymbol()
        count = 0
        for i, a in enumerate(mol.GetAtoms()):
            if a.GetSymbol() == symbol:
                count += 1
            if i == target_idx:
                return f"{symbol}{count}"
        return None

    def find_terminal_atom(self, smiles, smarts_pattern):
        """SMARTSパターンを用いて結合原子を自動特定"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: raise ValueError("Invalid SMILES.")
        pattern = Chem.MolFromSmarts(smarts_pattern)
        matches = mol.GetSubstructMatches(pattern)
        if not matches:
            raise ValueError(f"Pattern '{smarts_pattern}' not found.")
        return self._get_af3_atom_name(mol, matches[0][0])

    def prepare_af3_input(self, protein_seq, taca_smiles, bond_res_idx, terminal_smarts, bond_atom_protein="NZ"):
        """AlphaFold 3用JSONを作成"""
        bond_atom_ligand = self.find_terminal_atom(taca_smiles, terminal_smarts)
        data = {
            "name": self.job_name,
            "modelSeeds": [1],
            "sequences": [
                {"protein": {"id": "A", "sequence": protein_seq}},
                {"ligand": {"id": "B", "smiles": taca_smiles}}
            ],
            "bondedAtomPairs": [{
                "at1": {"resChainId": "A", "resIdx": bond_res_idx, "atomName": bond_atom_protein},
                "at2": {"resChainId": "B", "resIdx": 1, "atomName": bond_atom_ligand}
            }]
        }
        output_path = f"{self.job_name}.json"
        with open(output_path, "w") as f:
            json.dump(data, f, indent=4)
        return output_path

    def analyze_interactions(self, cif_path, distance_cutoff=5.0):
        """抗原内の接触残基を特定"""
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("model", cif_path)[0]
        protein_atoms = [a for a in structure['A'].get_atoms()]
        sugar_atoms = [a for a in structure.get_atoms() if "H_" in a.get_parent().get_full_id()[3][0]]
        ns = NeighborSearch(protein_atoms)
        contacts = {ns.search(s.coord, distance_cutoff, level="R")[0] for s in sugar_atoms if ns.search(s.coord, distance_cutoff, level="R")}
        return contacts

    def _calculate_sasa(self, cif_path):
        """FreeSASAを用いて糖鎖の露出面積を計算"""
        structure = freesasa.Structure(cif_path)
        result = freesasa.calc(structure)
        s_glyc = freesasa.selectObjects(["glyc, not chain A"], structure, result)
        return {"glycan_sasa": s_glyc["glyc"]}


class AntibodyDockingWorkflow:
    """抗体－抗原ドッキングとパラトープ解析を担当するクラス"""
    def __init__(self, job_name, h_chain="H", l_chain="L"):
        self.job_name = job_name
        self.h_chain = h_chain
        self.l_chain = l_chain

    def _get_cdr_labels(self, sequence, scheme="imgt"):
        """ANARCIによるCDRマッピング"""
        if anarci is None: return {}
        try:
            results = anarci([("seq", sequence)], scheme=scheme)
            numbering = results[0][0]
            if numbering is None: return {}
            labels = {}
            for pos_data, aa in numbering[0]:
                num = pos_data[0]
                label = "FW"
                if 27 <= num <= 38: label = "CDR1"
                elif 56 <= num <= 65: label = "CDR2"
                elif 105 <= num <= 117: label = "CDR3"
                labels[num] = label
            return labels
        except: return {}

    def analyze_paratope(self, cif_file, distance_cutoff=4.5):
        """パラトープ抽出とCDR割り当て"""
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("complex", cif_file)[0]
        antigen_atoms = [a for c in structure if c.id not in [self.h_chain, self.l_chain] for a in c.get_atoms()]
        ns = NeighborSearch([a for c in [self.h_chain, self.l_chain] if c in structure for a in structure[c].get_atoms()])
        
        contact_res = set()
        for a_atom in antigen_atoms:
            for res in ns.search(a_atom.coord, distance_cutoff, level="R"):
                contact_res.add(res)

        report = []
        for chain_id in [self.h_chain, self.l_chain]:
            if chain_id not in structure: continue
            seq = "".join([seq1(r.get_resname()) for r in structure[chain_id] if "CA" in r])
            cdr_map = self._get_cdr_labels(seq)
            for res in sorted(list(contact_res), key=lambda x: x.id[1]):
                if res.get_parent().id == chain_id:
                    report.append({
                        "Chain": "Heavy" if chain_id == self.h_chain else "Light",
                        "ResNum": res.id[1],
                        "ResName": res.get_resname(),
                        "CDR_Region": cdr_map.get(res.id[1], "FW")
                    })
        return pd.DataFrame(report)


class LightweightHotSpotAnalyzer:
    """接触密度に基づくHot Spot予測クラス"""
    def __init__(self, file_handle, h_chain="H", l_chain="L", ant_chains=["A", "B"]):
        # Streamlitから渡されるfile_handle(io.BytesIO)に対応
        content = file_handle.read().decode('utf-8')
        import io
        fh = io.StringIO(content)
        parser = MMCIFParser(QUIET=True)
        self.structure = parser.get_structure("complex", fh)[0]
        self.h_chain, self.l_chain, self.ant_chains = h_chain, l_chain, ant_chains

    def run_contact_density_scan(self, distance_cutoff=4.0):
        antigen_atoms = [a for cid in self.ant_chains if cid in self.structure for a in self.structure[cid].get_atoms()]
        ns = NeighborSearch(antigen_atoms)
        results = []
        for cid in [self.h_chain, self.l_chain]:
            if cid not in self.structure: continue
            for res in self.structure[cid]:
                if "CA" not in res: continue
                contact_count = sum(len(ns.search(a.coord, distance_cutoff)) for a in res.get_atoms())
                if contact_count > 0:
                    results.append({
                        "Chain": cid,
                        "Residue": f"{res.get_resname()}{res.id[1]}",
                        "HotSpot_Score": contact_count
                    })
        return pd.DataFrame(results).sort_values("HotSpot_Score", ascending=False) if results else pd.DataFrame()