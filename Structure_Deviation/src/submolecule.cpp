#include "submolecule.h"

atom::atom(string type, long int serial, string atom_name, string alt_ind, string res_name, string chain_name, long int res_no, double x, double y, double z, double occu, double t_fact, double radius, string atom_symbol, double partial_charge): Sphere(x, y, z, radius){
	_type = type;
	_serial_no = serial;
	_atom_name=atom_name;
	_alt_conf_indicator=alt_ind;
	_residue_name=res_name;
	_chain_name=chain_name;
	_residue_no=res_no;
	_occupancy = occu;
	_temp_factor = t_fact;
	_atom_symbol=atom_symbol;
	_partial_charge=partial_charge;

	assign_molecule_and_residue_type();
}

void atom::fprint(FILE* fp){
	fprintf(fp, "\n%s %ld %s %s %s %s %ld %lf %lf %lf %lf %lf %s %lf", _type.c_str(), _serial_no, _atom_name.c_str(), _alt_conf_indicator.c_str(), _residue_name.c_str(), _chain_name.c_str(), _residue_no, x, y, z, _occupancy, _temp_factor, _atom_symbol.c_str(), _partial_charge);
}

void atom::set_partial_charge(double pc){
	_partial_charge=pc;
}

void atom::set_atom_type(string at){
	_type=at;
}

void atom::set_serial_no(long int sno){
	_serial_no=sno;
}

void atom::set_atom_name(string an){
	_atom_name=an;
}

void atom::set_alt_conf_indicator(string ai){
	_alt_conf_indicator=ai;
}

void atom::set_residue_name(string rn){
	_residue_name=rn;
}

void atom::set_chain_name(string cn){
	_chain_name=cn;
}

void atom::set_residue_no(long int rno){
	_residue_no=rno;
}

void atom::set_occupancy(double oc){
	_occupancy=oc;
}

void atom::set_temp_factor(double tf){
	_temp_factor=tf;
}

void atom::set_atom_symbol(string as){
	_atom_symbol=as;
}

void atom::set_molecule_type(string mt){
	_molecule_type=mt;
}

void atom::set_residue_type(string rt){
	_residue_type=rt;
}

void atom::set_X(double x){
	_x_coord=x;
}

void atom::set_Y(double y){
	_y_coord=y;
}

void atom::set_Z(double z){
	_z_coord=z;
}

double atom::get_partial_charge(){
	return _partial_charge;
}

string atom::get_atom_type(){
	return _type;
}

long int atom::get_serial_no(){
	return _serial_no;
}

string atom::get_atom_name(){
	return _atom_name;
}

string atom::get_alt_conf_indicator(){
	return _alt_conf_indicator;
}

string atom::get_residue_name(){
	return _residue_name;
}

string atom::get_chain_name(){
	return _chain_name;
}

long int atom::get_residue_no(){
	return _residue_no;
}

double atom::get_occupancy(){
	return _occupancy;
}

double atom::get_temp_factor(){
	return _temp_factor;
}

string atom::get_atom_symbol(){
	return _atom_symbol;
}

string atom::get_molecule_type(){
	return _molecule_type;
}

string atom::get_residue_type(){
	return _residue_type;
}

void atom::assign_molecule_and_residue_type(){
	vector<string> standard_amino_acids {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
	vector<string> standard_dna {"DA", "DT", "DG", "DC", "DI", "DU"};
	vector<string> standard_rna {"A", "U", "G", "C", "I"};
	vector<string> non_standard_amino_acids {"TAV", "PHD", "OXX", "OHS", "KKD", "DMK", "D2T", "BHD", "BH2", "BFD", "ASL", "ASB", "AEI", "999", "0TD", "0AK", "0A0", "MGG", "ARO", "AGM", "2OR", "SD4", "MEN", "DMH", "AS7", "AHB", "AFA", "9DN", "823", "02L", "PVH", "2HF", "FZN", "LYX", "6CL", "RE3", "RE0", "R4K", "QMM", "NLQ", "MGN", "MFN", "MEQ", "LMQ", "HOO", "GNC", "GHG", "TS9", "QIL", "MEG", "LME", "IML", "ILX", "ILM", "IIL", "I2M", "G8M", "BIU", "S", "12L", "0E5", "05N", "037", "ZIQ", "TTQ", "TRX", "TRO", "TQI", "TOQ", "TNQ", "TGH", "TCR", "PAT", "O7D", "HTR", "HRP", "FTR", "FT6", "F7W", "E9M", "CTE", "BTR", "6CW", "5CW", "4PQ", "4IN", "4HT", "4FW", "1TQ", "0UO", "0AF", "WVL", "TBG", "PG1", "O7G", "MVA", "M2S", "LVN", "LE1", "IGL", "HVA", "FVA", "4SJ", "34E", "12L", "0E5", "0AB", "05N", "037", "033", "YOF", "U2X", "TYS", "TYQ", "TYN", "TYJ", "TYI", "TY9", "TY8", "TY5", "TY2", "TY1", "T0I", "PTR", "PTM", "PTH", "P3Q", "P2Q", "OMY", "OMX", "NIY", "NBQ", "MBQ", "IYR", "FY3", "FY2", "FLT", "F2Y", "DI7", "DBY", "DAH", "BYR", "AZY", "51T", "4LZ", "3YM", "3NF", "3MY", "3CT", "2TY", "2R3", "2LT", "1OP", "0WZ", "0PR", "0EA", "0A1", "ZU0", "Z3E", "YTH", "VAH", "TPO", "TMD", "THC", "TH6", "TH5", "T8L", "R4K", "QNY", "OTH", "OLT", "O7A", "NZC", "HY3", "HLU", "HL2", "H5M", "G3M", "EUP", "CTH", "BMT", "ALO", "8RE", "8JB", "6BR", "5CS", "4D4", "3PX", "2QZ", "28X", "26B", "0AH", "XDT", "WRP", "WPA", "UF0", "UDS", "TIS", "SXE", "SVZ", "SVY", "SVX", "SVW", "SVV", "SVA", "SUN", "SRZ", "SNM", "SGB", "SEP", "SEN", "SEM", "SEE", "SEB", "SDP", "SBL", "SBG", "SAC", "S1H", "S12", "RZ4", "RVX", "QDS", "P5U", "OSE", "OMH", "OLZ", "OAS", "NC1", "N10", "MIS", "MIR", "LPS", "LLY", "K5L", "J8W", "H14", "GVL", "GFT", "FGP", "FGL", "DDZ", "DBS", "BXT", "BG1", "BB8", "AZS", "ALS", "8SP", "73C", "5R5", "5JP", "57M", "54C", "4OZ", "4HJ", "4OV", "4HH", "432", "42Y", "3PX", "2ZC", "2RX", "1X6", "143", "DEA", "ZYK", "ZYJ", "YPR", "XPR", "VH0", "TPJ", "RT0", "PXU", "PRJ", "PR7", "PR4", "POM", "PLJ", "PH6", "PCA", "N80", "N7P", "MP8", "LWY", "JKH", "HZP", "HYP", "HY3", "H5M", "FPK", "FP9", "E0Y", "DYJ", "DPL", "6Y9", "4N9", "4N8", "4N7", "4L0", "4KY", "4FB", "45F", "3WX", "3PX", "3BY", "2P0", "12Y", "12X", "12L", "11Q", "0Y8", "0LF", "05N", "04U", "037", "SME", "SIB", "SAH", "OMT", "MT2", "MSL", "MME", "MHO", "ME0", "HYI", "HTI", "H1D", "FME", "ESC", "EPM", "CXM", "AME", "4MM", "4CY", "2FM", "XYC", "WLU", "TYT", "PAQ", "MNL", "MLL", "MLE", "MHL", "LEH", "LEF", "LED", "LAY", "HLU", "HL5", "HL2", "H7V", "G5G", "FLE", "F7Q", "CGU", "BL2", "ALC", "66D", "5GM", "2ML", "1L1", "11W", "0AG", "ZCL", "WPA", "WFP", "U3X", "TFQ", "TEF", "T11", "SMF", "QPH", "PPN", "PM3", "PHI", "PFF", "PF5", "PBF", "OZW", "NAL", "NA8", "N0A", "MTY", "MHU", "MEA", "IAM", "HOX", "H14", "FTY", "FCL", "FC0", "F2F", "DJD", "BNN", "BIF", "BB8", "ALN", "9NF", "7XC", "7N8", "6CV", "5CR", "4SJ", "4PH", "4OU", "4II", "4CF", "4BF", "4AF", "41H", "3CF", "33S", "2GX", "200", "1PA", "0BN", "0A9", "X2W", "SHR", "KGC", "GME", "GHW", "GHC", "G01", "EME", "CGA", "9NE", "3O3", "3GL", "11W", "ZZD", "ZBZ", "YCM", "XCN", "TSY", "TQZ", "TNB", "SNC", "SMC", "SLZ", "SIC", "SHC", "SCY", "SCS", "SCH", "S2C", "R1A", "QPA", "QCS", "PYX", "PRS", "PEC", "PBB", "P9S", "P1L", "OCY", "OCS", "NYS", "NYB", "NPH", "ML3", "MD5", "MD3", "MCS", "M2L", "M0H", "KPY", "KOR", "KNB", "K5H", "K1R", "JJL", "JJK", "JJJ", "J3D", "ICY", "HTI", "HNC", "GT9", "FOE", "FFM", "EJA", "EFC", "ECX", "DYS", "DC2", "CZZ", "CZ2", "CYW", "CYR", "CYQ", "CYG", "CYD", "CY4", "CY1", "CY0", "CSZ", "CSX", "CSU", "CSS", "CSR", "CSP", "CSO", "CSJ", "CSD", "CSB", "CSA", "CS4", "CS3", "CS1", "CMT", "CML", "CMH", "CME", "CGV", "CG6", "CCS", "CAS", "CAF", "C6C", "C5C", "C4R", "C3Y", "C1T", "C1S", "BUC", "BCS", "AGT", "8JB", "85F", "6WK", "6V1", "6M6", "60F", "5OW", "5CS", "4J4", "4GJ", "3X9", "31Q", "30V", "2XA", "2MT", "2CO", "143", "0QL", "0CS", "0A8", "07O", "03Y", "Z70", "Z01", "YPZ", "XX1", "XW1", "XOK", "VR0", "VB1", "UMA", "UM1", "TYY", "TTS", "TRQ", "TRN", "TRF", "TQQ", "TPQ", "TOX", "TLY", "TIH", "TCQ", "T9E", "T01", "SLL", "SKH", "SDB", "S2P", "RVJ", "RPI", "Q3P", "PYA", "PSH", "PRK", "POK", "PE1", "OYL", "ORQ", "ORN", "ONH", "OLD", "OHI", "OBS", "NWD", "NVA", "NPI", "NOT", "NNH", "NMM", "NLO", "NLG", "NLE", "NLB", "NEP", "NCB", "N9P", "N65", "MYN", "MYK", "MSO", "MSE", "MMO", "MLZ", "MLY", "MK8", "MHS", "MAA", "MA", "M3R", "M3L", "M2S", "02K", "02O", "02Y", "0AR", "0FL", "0X9", "AZH", "AYA", "AVJ", "API", "AN8", "AN6", "ALY", "AIB", "AHP", "AHO", "AGQ", "ABA", "AA4", "A8E", "A5N", "9WV", "9U0", "9TX", "9TU", "9TR", "9NV", "9NR", "9KP", "9E7", "8WY", "8LJ", "8AY", "7O5", "7JA", "74P", "73P", "73O", "73N", "6HN", "6GL", "6G4", "6DN", "5T3", "5OH", "5MW", "5GG", "5CT", "5AB", "56A", "4WQ", "4U7", "4J5", "4HL", "4DP", "4AW", "4AK", "41Q", "3ZH", "3WS", "3QN", "3AH", "33W", "2SO", "2RA", "2MR", "2KP", "2KK", "2AG", "23P", "1TY", "1AC", "B3U", "BP5", "BTK", "BWV", "C1J", "C1X", "C22", "C4G", "C67", "C6D", "CIR", "CLG", "CLH", "CWD", "CYJ", "CZS", "DA2", "DAB", "DBZ", "DDE", "DIR", "DLS", "DM0", "DNP", "DNS", "DNW", "DON", "DPP", "DPQ", "E9C", "E9V", "ELY", "ESB", "EXA", "EXY", "FAK", "FB5", "FB6", "FF9", "FH7", "FHL", "FHO", "FIO", "FL6", "FLA", "FQA", "GLJ", "GPL", "HAR", "HCM", "HHK", "HIC", "HIP", "HIQ", "HIX", "HLY", "HPE", "HQA", "HRG", "HS8", "HSE", "HSK", "HSL", "I58", "IEL", "ILY", "IOR", "IZO", "J9Y", "JLP", "KCR", "KCX", "KEO", "KFP", "KHB", "KPF", "KPI", "KST", "KYN", "KYQ", "L5P", "LA2", "LAL", "LBY", "LCK", "LDH", "LET", "LGY", "LLO", "LLP", "LLZ", "LMF", "LP6", "LPG", "LRK", "LSO", "LVG", "LYF", "LYO", "LYP", "LYR", "LYU", "LYZ", "XYG", "GMO", "NRP", "NRQ", "NYG", "QFG", "QLG", "RC7", "SWG", "X9Q", "5SQ", "BJO", "C99", "CFY", "CH6", "CH7", "CLV", "CR2", "CR7", "CRQ", "CRU", "CRX", "DYG", "EYG", "GYC", "GYS", "KZ1", "KZ4", "KZ7", "KZG", "KZY"};
	vector<string> non_standard_dna {"1AP", "2BU", "6HB", "A3A", "A40", "A5L", "A47", "AD2", "ABS", "ABR", "AF2", "E", "MA7", "R", "TCY", "XUA", "Y", "0AP", "1CC", "1FC", "47C", "4U3", "4PC", "5HC", "5FC", "5CM", "5SE", "5PC", "94O", "B7C", "C38", "C37", "C34", "C46", "C45", "CAR", "CBR", "CFZ", "CFL", "CSL", "D00", "D4B", "DNR", "EXC", "GCK", "ME6", "NCU", "TCJ", "TC1", "XCY", "0AD", "63H", "63G", "7BG", "8PY", "8MG", "8FG", "AFG", "BGM", "DGP", "DGI", "DG8", "DFG", "FMG", "F4Q", "G47", "G49", "GFL", "GF2", "GDR", "GSS", "GSR", "HN0", "M1G", "TTP", "TDP", "DGT", "8OG", "ZDU", "UPE", "UFR", "UBI", "TTM", "TLC", "TLN", "TFE", "TED", "TAF", "T41", "T39", "T38", "SMT", "PDU", "P2T", "NTT", "NMT", "NMS", "LST", "LSH", "JDT", "GMU", "EIT", "EAN", "DUZ", "DRT", "BVP", "BOE", "ATL", "8DT", "77Y", "5HU", "2ST", "2OT", "2NT", "2GT", "2BT", "2AT", "18Q", "XUG", "RDG", "PGN", "P", "MRG", "MG1", "LDG", "LCG"};
	vector<string> non_standard_rna {"YYG", "12A", "2OM", "AMO", "00A", "5FA", "6MZ", "6MA", "6IA", "A2L", "A23", "A2M", "A3P", "A44", "A9Z", "A5O", "ADP", "AET", "AP7", "AMD", "AVC", "GOM", "LCA", "MA6", "MIA", "N6G", "N79", "O2Z", "RIA", "T6A", "TXP", "TXD", "V3L", "3AU", "5FU", "5BU", "9QV", "75B", "FNU", "OMU", "UBD", "UAR", "U36", "U31", "UR3", "UPV", "IMP", "10C", "4OC", "5IC", "5HM", "92F", "A5M", "B8T", "C2L", "C31", "C43", "C5L", "CCC", "CBV", "CSF", "LHH", "LV2", "M5M", "N5M", "RSQ", "RPC", "PMT", "102", "1MG", "2EG", "23G", "CG1", "G1G", "G2L", "G3A", "G48", "GAO", "G7M", "GRB", "KAG", "GTP", "GDP", "UTP", "UDP", "U2L", "IU", "5GP", "GMP", "PSU", "1MA", "2MG", "5MC", "2MU", "5MU", "7MG", "H2U", "M2G", "O2G", "OMC", "OMG", "YG", "MNU", "CNU", "85Y", "3ME", "125", "126", "127", "TPG", "PGP", "MGV", "MGQ"};

	for (auto res: standard_amino_acids){
		if (_residue_name==res){
			_molecule_type="amino_acid";
			_residue_type="standard";
			return;
		}
	}

	for (auto res: non_standard_amino_acids){
		if (_residue_name==res){
			_molecule_type="amino_acid";
			_residue_type="non_stand";
			return;
		}
	}

	for (auto res: standard_dna){
		if (_residue_name==res){
			_molecule_type="dna_molecule";
			_residue_type="standard";
			return;
		}
	}

	for (auto res: non_standard_dna){
		if (_residue_name==res){
			_molecule_type="dna_molecule";
			_residue_type="non_stand";
			return;
		}
	}

	for (auto res: standard_rna){
		if (_residue_name==res){
			if (_type=="ATOM"){
				_molecule_type="rna_molecule";
				_residue_type="standard";
				return;
			}
		}
	}

	for (auto res: non_standard_rna){
		if (_residue_name==res){
			_molecule_type="rna_molecule";
			_residue_type="non_stand";
			return;
		}
	}

	if (_residue_name!="HOH"){
		_molecule_type="other_biom";
		_residue_type="non_stand";
	}

	if (_residue_name=="HOH"){
		_molecule_type="water_molecule";
		_residue_type="non_stand";
	}
}

residue::residue(vector<atom> atm_gr){
	_atom_group=atm_gr;
}

int residue::get_no_of_atoms(){
	return (_atom_group.size());
}

vector<atom> residue::get_atoms(){
	return _atom_group;
}

bond::bond(atom atm1, atom atm2, string bond_type, bool v1, bool v2){
	_atom1=atm1;
	_atom2=atm2;
	_bond_type=bond_type;
	_visited_atom1=v1;
	_visited_atom2=v2;
	compute_bond_strength();
}

atom bond::get_bonded_atom1(){
	return _atom1;
}

atom bond::get_bonded_atom2(){
	return _atom2;
}

string bond::get_bond_type(){
	return _bond_type;
}

void bond::set_visiting_status_atom1(bool v){
	_visited_atom1=v;
}

bool bond::get_visiting_status_atom1(){
	return _visited_atom1;
}

void bond::set_visiting_status_atom2(bool v){
	_visited_atom2=v;
}

bool bond::get_visiting_status_atom2(){
	return _visited_atom2;
}

double bond::get_bond_length(){
	return _bond_length;
}

double bond::get_bond_strength(){
	return _bond_strength;
}

void bond::compute_bond_strength(){
	double charge1=_atom1.get_partial_charge();
	double charge2=_atom2.get_partial_charge();
	double dist=_atom1.dist_sqr(&_atom2);
	_bond_strength=((charge1*charge2)/dist);
}
