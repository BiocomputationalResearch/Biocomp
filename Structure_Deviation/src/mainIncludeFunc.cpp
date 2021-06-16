#include "mainIncludeFunc.h"
#include "globalMain.h"
#include "mainIncludeFiles.h"
#include <openbabel/obconversion.h>

bool testNumber(char *str){
	for(int i = 0; i < strlen(str); i++){
		if (str[i] < 48 || str[i] > 57){
			return false;
		}
	}

	return true;
}

string convert_to_mol2(string prot_name, string accn, string ext){
	ifstream ifs(prot_name);

	string output_mol_file_name=accn+".mol2";
	ofstream ofs(output_mol_file_name);

	OpenBabel::OBConversion conv(&ifs, &ofs);
	if(!conv.SetInAndOutFormats(ext.substr(1, ext.size()).c_str(),"MOL2"))
	{
		cout << "Formats not available\n";
		output_mol_file_name="NOT_CONVERTED";
	}else{
		int n = conv.Convert();
		cout << n << " molecules converted\n";
	}

	return output_mol_file_name;
}

void form_tripl_sites(vector<metal_binding_site> tripl, string metal_name, string accn){
	if (metal_name=="NA"){
		no_of_compared_tripl_na+=tripl.size();
		for (auto ts: tripl){
			triplStructNa tripl_site=triplStructNa(ts, first_time_tripl_na, accn);
			if (first_time_tripl_na==true){
				first_time_tripl_na=false;
			}

			_tripl_na.push_back(tripl_site.form_tripl_site(no_of_compared_tripl_na));
		}
	}else if (metal_name=="CA"){
		no_of_compared_tripl_ca+=tripl.size();
		for (auto ts: tripl){
			triplStructCa tripl_site=triplStructCa(ts, first_time_tripl_ca, accn);
			if (first_time_tripl_ca==true){
				first_time_tripl_ca=false;
			}

			_tripl_ca.push_back(tripl_site.form_tripl_site(no_of_compared_tripl_ca));
		}
	}else if (metal_name=="K"){
		no_of_compared_tripl_k+=tripl.size();
		for (auto ts: tripl){
			triplStructK tripl_site=triplStructK(ts, first_time_tripl_k, accn);
			if (first_time_tripl_k==true){
				first_time_tripl_k=false;
			}

			_tripl_k.push_back(tripl_site.form_tripl_site(no_of_compared_tripl_k));
		}
	}else if (metal_name=="MG"){
		no_of_compared_tripl_mg+=tripl.size();
		for (auto ts: tripl){
			triplStructMg tripl_site=triplStructMg(ts, first_time_tripl_mg, accn);
			if (first_time_tripl_mg==true){
				first_time_tripl_mg=false;
			}

			_tripl_mg.push_back(tripl_site.form_tripl_site(no_of_compared_tripl_mg));
		}
	}else if (metal_name=="MN"){
		no_of_compared_tripl_mn+=tripl.size();
		for (auto ts: tripl){
			triplStructMn tripl_site=triplStructMn(ts, first_time_tripl_mn, accn);
			if (first_time_tripl_mn==true){
				first_time_tripl_mn=false;
			}

			_tripl_mn.push_back(tripl_site.form_tripl_site(no_of_compared_tripl_mn));
		}
	}else if (metal_name=="FE"){
		no_of_compared_tripl_fe+=tripl.size();
		for (auto ts: tripl){
			triplStructFe tripl_site=triplStructFe(ts, first_time_tripl_fe, accn);
			if (first_time_tripl_fe==true){
				first_time_tripl_fe=false;
			}

			_tripl_fe.push_back(tripl_site.form_tripl_site(no_of_compared_tripl_fe));
		}
	}else if (metal_name=="CU"){
		no_of_compared_tripl_cu+=tripl.size();
		for (auto ts: tripl){
			triplStructCu tripl_site=triplStructCu(ts, first_time_tripl_cu, accn);
			if (first_time_tripl_cu==true){
				first_time_tripl_cu=false;
			}

			_tripl_cu.push_back(tripl_site.form_tripl_site(no_of_compared_tripl_cu));
		}
	}else if (metal_name=="ZN"){
		no_of_compared_tripl_zn+=tripl.size();
		for (auto ts: tripl){
			triplStructZn tripl_site=triplStructZn(ts, first_time_tripl_zn, accn);
			if (first_time_tripl_zn==true){
				first_time_tripl_zn=false;
			}

			_tripl_zn.push_back(tripl_site.form_tripl_site(no_of_compared_tripl_zn));
		}
	}else if (metal_name=="CO"){
		no_of_compared_tripl_co+=tripl.size();
		for (auto ts: tripl){
			triplStructCo tripl_site=triplStructCo(ts, first_time_tripl_co, accn);
			if (first_time_tripl_co==true){
				first_time_tripl_co=false;
			}

			_tripl_co.push_back(tripl_site.form_tripl_site(no_of_compared_tripl_co));
		}
	}else if (metal_name=="PB"){
		no_of_compared_tripl_pb+=tripl.size();
		for (auto ts: tripl){
			triplStructPb tripl_site=triplStructPb(ts, first_time_tripl_pb, accn);
			if (first_time_tripl_pb==true){
				first_time_tripl_pb=false;
			}

			_tripl_pb.push_back(tripl_site.form_tripl_site(no_of_compared_tripl_pb));
		}
	}else if (metal_name=="HG"){
		no_of_compared_tripl_hg+=tripl.size();
		for (auto ts: tripl){
			triplStructHg tripl_site=triplStructHg(ts, first_time_tripl_hg, accn);
			if (first_time_tripl_hg==true){
				first_time_tripl_hg=false;
			}

			_tripl_hg.push_back(tripl_site.form_tripl_site(no_of_compared_tripl_hg));
		}
	}else if (metal_name=="AS"){
		no_of_compared_tripl_as+=tripl.size();
		for (auto ts: tripl){
			triplStructAs tripl_site=triplStructAs(ts, first_time_tripl_as, accn);
			if (first_time_tripl_as==true){
				first_time_tripl_as=false;
			}

			_tripl_as.push_back(tripl_site.form_tripl_site(no_of_compared_tripl_as));
		}
	}else if (metal_name=="CD"){
		no_of_compared_tripl_cd+=tripl.size();
		for (auto ts: tripl){
			triplStructCd tripl_site=triplStructCd(ts, first_time_tripl_cd, accn);
			if (first_time_tripl_cd==true){
				first_time_tripl_cd=false;
			}

			_tripl_cd.push_back(tripl_site.form_tripl_site(no_of_compared_tripl_cd));
		}
	}else if (metal_name=="NI"){
		no_of_compared_tripl_ni+=tripl.size();
		for (auto ts: tripl){
			triplStructNi tripl_site=triplStructNi(ts, first_time_tripl_ni, accn);
			if (first_time_tripl_ni==true){
				first_time_tripl_ni=false;
			}

			_tripl_ni.push_back(tripl_site.form_tripl_site(no_of_compared_tripl_ni));
		}
	}
}

void form_tetrapl_sites(vector<metal_binding_site> tetrapl, string metal_name, string accn){
	if (metal_name=="MG"){
		no_of_compared_tetrapl_mg+=tetrapl.size();
		for (auto ts: tetrapl){
			tetraplStructMg tetrapl_site=tetraplStructMg(ts, first_time_tetrapl_mg, accn);
			if (first_time_tetrapl_mg==true){
				first_time_tetrapl_mg=false;
			}

			_tetrapl_mg.push_back(tetrapl_site.form_tetrapl_site(no_of_compared_tetrapl_mg));
		}
	}else if (metal_name=="FE"){
		no_of_compared_tetrapl_fe+=tetrapl.size();
		for (auto ts: tetrapl){
			tetraplStructFe tetrapl_site=tetraplStructFe(ts, first_time_tetrapl_fe, accn);
			if (first_time_tetrapl_fe==true){
				first_time_tetrapl_fe=false;
			}

			_tetrapl_fe.push_back(tetrapl_site.form_tetrapl_site(no_of_compared_tetrapl_fe));
		}
	}else if (metal_name=="CD"){
		no_of_compared_tetrapl_cd+=tetrapl.size();
		for (auto ts: tetrapl){
			tetraplStructCd tetrapl_site=tetraplStructCd(ts, first_time_tetrapl_cd, accn);
			if (first_time_tetrapl_cd==true){
				first_time_tetrapl_cd=false;
			}

			_tetrapl_cd.push_back(tetrapl_site.form_tetrapl_site(no_of_compared_tetrapl_cd));
		}
	}else if (metal_name=="NI"){
		no_of_compared_tetrapl_ni+=tetrapl.size();
		for (auto ts: tetrapl){
			tetraplStructNi tetrapl_site=tetraplStructNi(ts, first_time_tetrapl_ni, accn);
			if (first_time_tetrapl_ni==true){
				first_time_tetrapl_ni=false;
			}

			_tetrapl_ni.push_back(tetrapl_site.form_tetrapl_site(no_of_compared_tetrapl_ni));
		}
	}
}

void form_tetra_sites(vector<metal_binding_site> tetra, string metal_name, string accn){
	if (metal_name=="NA"){
		no_of_compared_tetra_na+=tetra.size();
		for (auto ts: tetra){
			tetraStructNa tetra_site=tetraStructNa(ts, first_time_tetra_na, accn);
			if (first_time_tetra_na==true){
				first_time_tetra_na=false;
			}

			_tetra_na.push_back(tetra_site.form_tetra_site(no_of_compared_tetra_na));
		}
	}else if (metal_name=="CA"){
		no_of_compared_tetra_ca+=tetra.size();
		for (auto ts: tetra){
			tetraStructCa tetra_site=tetraStructCa(ts, first_time_tetra_ca, accn);
			if (first_time_tetra_ca==true){
				first_time_tetra_ca=false;
			}

			_tetra_ca.push_back(tetra_site.form_tetra_site(no_of_compared_tetra_ca));
		}
	}else if (metal_name=="K"){
		no_of_compared_tetra_k+=tetra.size();
		for (auto ts: tetra){
			tetraStructK tetra_site=tetraStructK(ts, first_time_tetra_k, accn);
			if (first_time_tetra_k==true){
				first_time_tetra_k=false;
			}

			_tetra_k.push_back(tetra_site.form_tetra_site(no_of_compared_tetra_k));
		}
	}else if (metal_name=="MG"){
		no_of_compared_tetra_mg+=tetra.size();
		for (auto ts: tetra){
			tetraStructMg tetra_site=tetraStructMg(ts, first_time_tetra_mg, accn);
			if (first_time_tetra_mg==true){
				first_time_tetra_mg=false;
			}

			_tetra_mg.push_back(tetra_site.form_tetra_site(no_of_compared_tetra_mg));
		}
	}else if (metal_name=="MN"){
		no_of_compared_tetra_mn+=tetra.size();
		for (auto ts: tetra){
			tetraStructMn tetra_site=tetraStructMn(ts, first_time_tetra_mn, accn);
			if (first_time_tetra_mn==true){
				first_time_tetra_mn=false;
			}

			_tetra_mn.push_back(tetra_site.form_tetra_site(no_of_compared_tetra_mn));
		}
	}else if (metal_name=="FE"){
		no_of_compared_tetra_fe+=tetra.size();
		for (auto ts: tetra){
			tetraStructFe tetra_site=tetraStructFe(ts, first_time_tetra_fe, accn);
			if (first_time_tetra_fe==true){
				first_time_tetra_fe=false;
			}

			_tetra_fe.push_back(tetra_site.form_tetra_site(no_of_compared_tetra_fe));
		}
	}else if (metal_name=="CU"){
		no_of_compared_tetra_cu+=tetra.size();
		for (auto ts: tetra){
			tetraStructCu tetra_site=tetraStructCu(ts, first_time_tetra_cu, accn);
			if (first_time_tetra_cu==true){
				first_time_tetra_cu=false;
			}

			_tetra_cu.push_back(tetra_site.form_tetra_site(no_of_compared_tetra_cu));
		}
	}else if (metal_name=="ZN"){
		no_of_compared_tetra_zn+=tetra.size();
		for (auto ts: tetra){
			tetraStructZn tetra_site=tetraStructZn(ts, first_time_tetra_zn, accn);
			if (first_time_tetra_zn==true){
				first_time_tetra_zn=false;
			}

			_tetra_zn.push_back(tetra_site.form_tetra_site(no_of_compared_tetra_zn));
		}
	}else if (metal_name=="CO"){
		no_of_compared_tetra_co+=tetra.size();
		for (auto ts: tetra){
			tetraStructCo tetra_site=tetraStructCo(ts, first_time_tetra_co, accn);
			if (first_time_tetra_co==true){
				first_time_tetra_co=false;
			}

			_tetra_co.push_back(tetra_site.form_tetra_site(no_of_compared_tetra_co));
		}
	}else if (metal_name=="PB"){
		no_of_compared_tetra_pb+=tetra.size();
		for (auto ts: tetra){
			tetraStructPb tetra_site=tetraStructPb(ts, first_time_tetra_pb, accn);
			if (first_time_tetra_pb==true){
				first_time_tetra_pb=false;
			}

			_tetra_pb.push_back(tetra_site.form_tetra_site(no_of_compared_tetra_pb));
		}
	}else if (metal_name=="HG"){
		no_of_compared_tetra_hg+=tetra.size();
		for (auto ts: tetra){
			tetraStructHg tetra_site=tetraStructHg(ts, first_time_tetra_hg, accn);
			if (first_time_tetra_hg==true){
				first_time_tetra_hg=false;
			}

			_tetra_hg.push_back(tetra_site.form_tetra_site(no_of_compared_tetra_hg));
		}
	}else if (metal_name=="AS"){
		no_of_compared_tetra_as+=tetra.size();
		for (auto ts: tetra){
			tetraStructAs tetra_site=tetraStructAs(ts, first_time_tetra_as, accn);
			if (first_time_tetra_as==true){
				first_time_tetra_as=false;
			}

			_tetra_as.push_back(tetra_site.form_tetra_site(no_of_compared_tetra_as));
		}
	}else if (metal_name=="CD"){
		no_of_compared_tetra_cd+=tetra.size();
		for (auto ts: tetra){
			tetraStructCd tetra_site=tetraStructCd(ts, first_time_tetra_cd, accn);
			if (first_time_tetra_cd==true){
				first_time_tetra_cd=false;
			}

			_tetra_cd.push_back(tetra_site.form_tetra_site(no_of_compared_tetra_cd));
		}
	}else if (metal_name=="NI"){
		no_of_compared_tetra_ni+=tetra.size();
		for (auto ts: tetra){
			tetraStructNi tetra_site=tetraStructNi(ts, first_time_tetra_ni, accn);
			if (first_time_tetra_ni==true){
				first_time_tetra_ni=false;
			}

			_tetra_ni.push_back(tetra_site.form_tetra_site(no_of_compared_tetra_ni));
		}
	}else if (metal_name=="CS"){
		no_of_compared_tetra_cs+=tetra.size();
		for (auto ts: tetra){
			tetraStructCs tetra_site=tetraStructCs(ts, first_time_tetra_cs, accn);
			if (first_time_tetra_cs==true){
				first_time_tetra_cs=false;
			}

			_tetra_cs.push_back(tetra_site.form_tetra_site(no_of_compared_tetra_cs));
		}
	}
}

void form_tetrapyr_sites(vector<metal_binding_site> tetrapyr, string metal_name, string accn){
	if (metal_name=="NA"){
		no_of_compared_tetrapyr_na+=tetrapyr.size();
		for (auto tp: tetrapyr){
			tetrapyrStructNa tetrapyr_site=tetrapyrStructNa(tp, first_time_tetrapyr_na, accn);
			if (first_time_tetrapyr_na==true){
				first_time_tetrapyr_na=false;
			}

			_tetrapyr_na.push_back(tetrapyr_site.form_tetrapyr_site(no_of_compared_tetrapyr_na));
		}
	}else if (metal_name=="CA"){
		no_of_compared_tetrapyr_ca+=tetrapyr.size();
		for (auto tp: tetrapyr){
			tetrapyrStructCa tetrapyr_site=tetrapyrStructCa(tp, first_time_tetrapyr_ca, accn);
			if (first_time_tetrapyr_ca==true){
				first_time_tetrapyr_ca=false;
			}

			_tetrapyr_ca.push_back(tetrapyr_site.form_tetrapyr_site(no_of_compared_tetrapyr_ca));
		}
	}else if (metal_name=="MG"){
		no_of_compared_tetrapyr_mg+=tetrapyr.size();
		for (auto tp: tetrapyr){
			tetrapyrStructMg tetrapyr_site=tetrapyrStructMg(tp, first_time_tetrapyr_mg, accn);
			if (first_time_tetrapyr_mg==true){
				first_time_tetrapyr_mg=false;
			}

			_tetrapyr_mg.push_back(tetrapyr_site.form_tetrapyr_site(no_of_compared_tetrapyr_mg));
		}
	}else if (metal_name=="CD"){
		no_of_compared_tetrapyr_cd+=tetrapyr.size();
		for (auto tp: tetrapyr){
			tetrapyrStructCd tetrapyr_site=tetrapyrStructCd(tp, first_time_tetrapyr_cd, accn);
			if (first_time_tetrapyr_cd==true){
				first_time_tetrapyr_cd=false;
			}

			_tetrapyr_cd.push_back(tetrapyr_site.form_tetrapyr_site(no_of_compared_tetrapyr_cd));
		}
	}else if (metal_name=="NI"){
		no_of_compared_tetrapyr_ni+=tetrapyr.size();
		for (auto tp: tetrapyr){
			tetrapyrStructNi tetrapyr_site=tetrapyrStructNi(tp, first_time_tetrapyr_ni, accn);
			if (first_time_tetrapyr_ni==true){
				first_time_tetrapyr_ni=false;
			}

			_tetrapyr_ni.push_back(tetrapyr_site.form_tetrapyr_site(no_of_compared_tetrapyr_ni));
		}
	}else if (metal_name=="ZN"){
		no_of_compared_tetrapyr_zn+=tetrapyr.size();
		for (auto tp: tetrapyr){
			tetrapyrStructZn tetrapyr_site=tetrapyrStructZn(tp, first_time_tetrapyr_zn, accn);
			if (first_time_tetrapyr_zn==true){
				first_time_tetrapyr_zn=false;
			}

			_tetrapyr_zn.push_back(tetrapyr_site.form_tetrapyr_site(no_of_compared_tetrapyr_zn));
		}
	}
}

void form_tribipyr_sites(vector<metal_binding_site> tribipyr, string metal_name, string accn){
	if (metal_name=="NA"){
		no_of_compared_tribipyr_na+=tribipyr.size();
		for (auto tbp: tribipyr){
			tribipyrStructNa tribipyr_site=tribipyrStructNa(tbp, first_time_tribipyr_na, accn);
			if (first_time_tribipyr_na==true){
				first_time_tribipyr_na=false;
			}

			_tribipyr_na.push_back(tribipyr_site.form_tribipyr_site(no_of_compared_tribipyr_na));
		}
	}else if (metal_name=="CA"){
		no_of_compared_tribipyr_ca+=tribipyr.size();
		for (auto tbp: tribipyr){
			tribipyrStructCa tribipyr_site=tribipyrStructCa(tbp, first_time_tribipyr_ca, accn);
			if (first_time_tribipyr_ca==true){
				first_time_tribipyr_ca=false;
			}

			_tribipyr_ca.push_back(tribipyr_site.form_tribipyr_site(no_of_compared_tribipyr_ca));
		}
	}else if (metal_name=="MG"){
		no_of_compared_tribipyr_mg+=tribipyr.size();
		for (auto tbp: tribipyr){
			tribipyrStructMg tribipyr_site=tribipyrStructMg(tbp, first_time_tribipyr_mg, accn);
			if (first_time_tribipyr_mg==true){
				first_time_tribipyr_mg=false;
			}

			_tribipyr_mg.push_back(tribipyr_site.form_tribipyr_site(no_of_compared_tribipyr_mg));
		}
	}else if (metal_name=="MN"){
		no_of_compared_tribipyr_mn+=tribipyr.size();
		for (auto tbp: tribipyr){
			tribipyrStructMn tribipyr_site=tribipyrStructMn(tbp, first_time_tribipyr_mn, accn);
			if (first_time_tribipyr_mn==true){
				first_time_tribipyr_mn=false;
			}

			_tribipyr_mn.push_back(tribipyr_site.form_tribipyr_site(no_of_compared_tribipyr_mn));
		}
	}else if (metal_name=="FE"){
		no_of_compared_tribipyr_fe+=tribipyr.size();
		for (auto tbp: tribipyr){
			tribipyrStructFe tribipyr_site=tribipyrStructFe(tbp, first_time_tribipyr_fe, accn);
			if (first_time_tribipyr_fe==true){
				first_time_tribipyr_fe=false;
			}

			_tribipyr_fe.push_back(tribipyr_site.form_tribipyr_site(no_of_compared_tribipyr_fe));
		}
	}else if (metal_name=="ZN"){
		no_of_compared_tribipyr_zn+=tribipyr.size();
		for (auto tbp: tribipyr){
			tribipyrStructZn tribipyr_site=tribipyrStructZn(tbp, first_time_tribipyr_zn, accn);
			if (first_time_tribipyr_zn==true){
				first_time_tribipyr_zn=false;
			}

			_tribipyr_zn.push_back(tribipyr_site.form_tribipyr_site(no_of_compared_tribipyr_zn));
		}
	}else if (metal_name=="CO"){
		no_of_compared_tribipyr_co+=tribipyr.size();
		for (auto tbp: tribipyr){
			tribipyrStructCo tribipyr_site=tribipyrStructCo(tbp, first_time_tribipyr_co, accn);
			if (first_time_tribipyr_co==true){
				first_time_tribipyr_co=false;
			}

			_tribipyr_co.push_back(tribipyr_site.form_tribipyr_site(no_of_compared_tribipyr_co));
		}
	}else if (metal_name=="HG"){
		no_of_compared_tribipyr_hg+=tribipyr.size();
		for (auto tbp: tribipyr){
			tribipyrStructHg tribipyr_site=tribipyrStructHg(tbp, first_time_tribipyr_hg, accn);
			if (first_time_tribipyr_hg==true){
				first_time_tribipyr_hg=false;
			}

			_tribipyr_hg.push_back(tribipyr_site.form_tribipyr_site(no_of_compared_tribipyr_hg));
		}
	}else if (metal_name=="CD"){
		no_of_compared_tribipyr_cd+=tribipyr.size();
		for (auto tbp: tribipyr){
			tribipyrStructCd tribipyr_site=tribipyrStructCd(tbp, first_time_tribipyr_cd, accn);
			if (first_time_tribipyr_cd==true){
				first_time_tribipyr_cd=false;
			}

			_tribipyr_cd.push_back(tribipyr_site.form_tribipyr_site(no_of_compared_tribipyr_cd));
		}
	}else if (metal_name=="NI"){
		no_of_compared_tribipyr_ni+=tribipyr.size();
		for (auto tbp: tribipyr){
			tribipyrStructNi tribipyr_site=tribipyrStructNi(tbp, first_time_tribipyr_ni, accn);
			if (first_time_tribipyr_ni==true){
				first_time_tribipyr_ni=false;
			}

			_tribipyr_ni.push_back(tribipyr_site.form_tribipyr_site(no_of_compared_tribipyr_ni));
		}
	}
}

void form_octa_sites(vector<metal_binding_site> octa, string metal_name, string accn){
	cout<<"\nForm Octahedron";
	if (metal_name=="NA"){
		no_of_compared_octa_na+=octa.size();
		for (auto oct: octa){
			octaStructNa octa_site=octaStructNa(oct, first_time_octa_na, accn);
			if (first_time_octa_na==true){
				first_time_octa_na=false;
			}

			_octa_na.push_back(octa_site.form_octa_site(no_of_compared_octa_na));
		}
	}else if (metal_name=="CA"){
		no_of_compared_octa_ca+=octa.size();
		for (auto oct: octa){
			octaStructCa octa_site=octaStructCa(oct, first_time_octa_ca, accn);
			if (first_time_octa_ca==true){
				first_time_octa_ca=false;
			}

			_octa_ca.push_back(octa_site.form_octa_site(no_of_compared_octa_ca));
		}
	}else if (metal_name=="MG"){
		no_of_compared_octa_mg+=octa.size();
		for (auto oct: octa){
			octaStructMg octa_site=octaStructMg(oct, first_time_octa_mg, accn);
			if (first_time_octa_mg==true){
				first_time_octa_mg=false;
			}

			_octa_mg.push_back(octa_site.form_octa_site(no_of_compared_octa_mg));
		}
	}else if (metal_name=="MN"){
		no_of_compared_octa_mn+=octa.size();
		for (auto oct: octa){
			octaStructMn octa_site=octaStructMn(oct, first_time_octa_mn, accn);
			if (first_time_octa_mn==true){
				first_time_octa_mn=false;
			}

			_octa_mn.push_back(octa_site.form_octa_site(no_of_compared_octa_mn));
		}
	}else if (metal_name=="FE"){
		no_of_compared_octa_fe+=octa.size();
		for (auto oct: octa){
			octaStructFe octa_site=octaStructFe(oct, first_time_octa_fe, accn);
			if (first_time_octa_fe==true){
				first_time_octa_fe=false;
			}

			_octa_fe.push_back(octa_site.form_octa_site(no_of_compared_octa_fe));
		}
	}else if (metal_name=="CD"){
		no_of_compared_octa_cd+=octa.size();
		for (auto oct: octa){
			octaStructCd octa_site=octaStructCd(oct, first_time_octa_cd, accn);
			if (first_time_octa_cd==true){
				first_time_octa_cd=false;
			}

			_octa_cd.push_back(octa_site.form_octa_site(no_of_compared_octa_cd));
		}
	}else if (metal_name=="NI"){
		no_of_compared_octa_ni+=octa.size();
		for (auto oct: octa){
			octaStructNi octa_site=octaStructNi(oct, first_time_octa_ni, accn);
			if (first_time_octa_ni==true){
				first_time_octa_ni=false;
			}

			_octa_ni.push_back(octa_site.form_octa_site(no_of_compared_octa_ni));
		}
	}else if (metal_name=="ZN"){
		no_of_compared_octa_zn+=octa.size();
		for (auto oct: octa){
			octaStructZn octa_site=octaStructZn(oct, first_time_octa_zn, accn);
			if (first_time_octa_zn==true){
				first_time_octa_zn=false;
			}

			_octa_zn.push_back(octa_site.form_octa_site(no_of_compared_octa_zn));
		}
	}else if (metal_name=="CO"){
		no_of_compared_octa_co+=octa.size();
		for (auto oct: octa){
			octaStructCo octa_site=octaStructCo(oct, first_time_octa_co, accn);
			if (first_time_octa_co==true){
				first_time_octa_co=false;
			}

			_octa_co.push_back(octa_site.form_octa_site(no_of_compared_octa_co));
		}
	}
}

void generate_summary_report_tripl(){
	if (no_of_compared_tripl_na!=0){
		triplStructNa::gen_summary_report_tripl_na(no_of_compared_tripl_na);
	}

	if (no_of_compared_tripl_ca!=0){
		triplStructCa::gen_summary_report_tripl_ca(no_of_compared_tripl_ca);
	}

	if (no_of_compared_tripl_k!=0){
		triplStructK::gen_summary_report_tripl_k(no_of_compared_tripl_k);
	}

	if (no_of_compared_tripl_mg!=0){
		triplStructMg::gen_summary_report_tripl_mg(no_of_compared_tripl_mg);
	}

	if (no_of_compared_tripl_mn!=0){
		triplStructMn::gen_summary_report_tripl_mn(no_of_compared_tripl_mn);
	}

	if (no_of_compared_tripl_fe!=0){
		triplStructFe::gen_summary_report_tripl_fe(no_of_compared_tripl_fe);
	}

	if (no_of_compared_tripl_cu!=0){
		triplStructCu::gen_summary_report_tripl_cu(no_of_compared_tripl_cu);
	}

	if (no_of_compared_tripl_zn!=0){
		triplStructZn::gen_summary_report_tripl_zn(no_of_compared_tripl_zn);
	}

	if (no_of_compared_tripl_co!=0){
		triplStructCo::gen_summary_report_tripl_co(no_of_compared_tripl_co);
	}

	if (no_of_compared_tripl_pb!=0){
		triplStructPb::gen_summary_report_tripl_pb(no_of_compared_tripl_pb);
	}

	if (no_of_compared_tripl_hg!=0){
		triplStructHg::gen_summary_report_tripl_hg(no_of_compared_tripl_hg);
	}

	if (no_of_compared_tripl_as!=0){
		triplStructAs::gen_summary_report_tripl_as(no_of_compared_tripl_as);
	}

	if (no_of_compared_tripl_cd!=0){
		triplStructCd::gen_summary_report_tripl_cd(no_of_compared_tripl_cd);
	}

	if (no_of_compared_tripl_ni!=0){
		triplStructNi::gen_summary_report_tripl_ni(no_of_compared_tripl_ni);
	}
}

void generate_summary_report_tetrapl(){
	if (no_of_compared_tetrapl_mg!=0){
		tetraplStructMg::gen_summary_report_tetrapl_mg(no_of_compared_tetrapl_mg);
	}

	if (no_of_compared_tetrapl_fe!=0){
		tetraplStructFe::gen_summary_report_tetrapl_fe(no_of_compared_tetrapl_fe);
	}

	if (no_of_compared_tetrapl_cd!=0){
		tetraplStructCd::gen_summary_report_tetrapl_cd(no_of_compared_tetrapl_cd);
	}

	if (no_of_compared_tetrapl_ni!=0){
		tetraplStructNi::gen_summary_report_tetrapl_ni(no_of_compared_tetrapl_ni);
	}
}

void generate_summary_report_tetra(){
	if (no_of_compared_tetra_na!=0){
		tetraStructNa::gen_summary_report_tetra_na(no_of_compared_tetra_na);
	}

	if (no_of_compared_tetra_ca!=0){
		tetraStructCa::gen_summary_report_tetra_ca(no_of_compared_tetra_ca);
	}

	if (no_of_compared_tetra_k!=0){
		tetraStructK::gen_summary_report_tetra_k(no_of_compared_tetra_k);
	}

	if (no_of_compared_tetra_mg!=0){
		tetraStructMg::gen_summary_report_tetra_mg(no_of_compared_tetra_mg);
	}

	if (no_of_compared_tetra_mn!=0){
		tetraStructMn::gen_summary_report_tetra_mn(no_of_compared_tetra_mn);
	}

	if (no_of_compared_tetra_fe!=0){
		tetraStructFe::gen_summary_report_tetra_fe(no_of_compared_tetra_fe);
	}

	if (no_of_compared_tetra_cu!=0){
		tetraStructCu::gen_summary_report_tetra_cu(no_of_compared_tetra_cu);
	}

	if (no_of_compared_tetra_zn!=0){
		tetraStructZn::gen_summary_report_tetra_zn(no_of_compared_tetra_zn);
	}

	if (no_of_compared_tetra_co!=0){
		tetraStructCo::gen_summary_report_tetra_co(no_of_compared_tetra_co);
	}

	if (no_of_compared_tetra_pb!=0){
		tetraStructPb::gen_summary_report_tetra_pb(no_of_compared_tetra_pb);
	}

	if (no_of_compared_tetra_hg!=0){
		tetraStructHg::gen_summary_report_tetra_hg(no_of_compared_tetra_hg);
	}

	if (no_of_compared_tetra_as!=0){
		tetraStructAs::gen_summary_report_tetra_as(no_of_compared_tetra_as);
	}

	if (no_of_compared_tetra_cd!=0){
		tetraStructCd::gen_summary_report_tetra_cd(no_of_compared_tetra_cd);
	}

	if (no_of_compared_tetra_ni!=0){
		tetraStructNi::gen_summary_report_tetra_ni(no_of_compared_tetra_ni);
	}

	if (no_of_compared_tetra_cs!=0){
		tetraStructCs::gen_summary_report_tetra_cs(no_of_compared_tetra_cs);
	}
}

void generate_summary_report_tetrapyr(){
	if (no_of_compared_tetrapyr_na!=0){
		tetrapyrStructNa::gen_summary_report_tetrapyr_na(no_of_compared_tetrapyr_na);
	}

	if (no_of_compared_tetrapyr_ca!=0){
		tetrapyrStructCa::gen_summary_report_tetrapyr_ca(no_of_compared_tetrapyr_ca);
	}

	if (no_of_compared_tetrapyr_mg!=0){
		tetrapyrStructMg::gen_summary_report_tetrapyr_mg(no_of_compared_tetrapyr_mg);
	}

	if (no_of_compared_tetrapyr_cd!=0){
		tetrapyrStructCd::gen_summary_report_tetrapyr_cd(no_of_compared_tetrapyr_cd);
	}

	if (no_of_compared_tetrapyr_ni!=0){
		tetrapyrStructNi::gen_summary_report_tetrapyr_ni(no_of_compared_tetrapyr_ni);
	}

	if (no_of_compared_tetrapyr_zn!=0){
		tetrapyrStructZn::gen_summary_report_tetrapyr_zn(no_of_compared_tetrapyr_zn);
	}
}

void generate_summary_report_octa(){
	if (no_of_compared_octa_na!=0){
		octaStructNa::gen_summary_report_octa_na(no_of_compared_octa_na);
	}

	if (no_of_compared_octa_ca!=0){
		octaStructCa::gen_summary_report_octa_ca(no_of_compared_octa_ca);
	}

	if (no_of_compared_octa_mg!=0){
		octaStructMg::gen_summary_report_octa_mg(no_of_compared_octa_mg);
	}

	if (no_of_compared_octa_mn!=0){
		octaStructMn::gen_summary_report_octa_mn(no_of_compared_octa_mn);
	}

	if (no_of_compared_octa_fe!=0){
		octaStructFe::gen_summary_report_octa_fe(no_of_compared_octa_fe);
	}

	if (no_of_compared_octa_cd!=0){
		octaStructCd::gen_summary_report_octa_cd(no_of_compared_octa_cd);
	}

	if (no_of_compared_octa_ni!=0){
		octaStructNi::gen_summary_report_octa_ni(no_of_compared_octa_ni);
	}

	if (no_of_compared_octa_zn!=0){
		octaStructZn::gen_summary_report_octa_zn(no_of_compared_octa_zn);
	}

	if (no_of_compared_octa_co!=0){
		octaStructCo::gen_summary_report_octa_co(no_of_compared_octa_co);
	}
}

void generate_summary_report_tribipyr(){
	if (no_of_compared_tribipyr_na!=0){
		tribipyrStructNa::gen_summary_report_tribipyr_na(no_of_compared_tribipyr_na);
	}

	if (no_of_compared_tribipyr_ca!=0){
		tribipyrStructCa::gen_summary_report_tribipyr_ca(no_of_compared_tribipyr_ca);
	}

	if (no_of_compared_tribipyr_mg!=0){
		tribipyrStructMg::gen_summary_report_tribipyr_mg(no_of_compared_tribipyr_mg);
	}

	if (no_of_compared_tribipyr_mn!=0){
		tribipyrStructMn::gen_summary_report_tribipyr_mn(no_of_compared_tribipyr_mn);
	}

	if (no_of_compared_tribipyr_fe!=0){
		tribipyrStructFe::gen_summary_report_tribipyr_fe(no_of_compared_tribipyr_fe);
	}

	if (no_of_compared_tribipyr_zn!=0){
		tribipyrStructZn::gen_summary_report_tribipyr_zn(no_of_compared_tribipyr_zn);
	}

	if (no_of_compared_tribipyr_co!=0){
		tribipyrStructCo::gen_summary_report_tribipyr_co(no_of_compared_tribipyr_co);
	}

	if (no_of_compared_tribipyr_hg!=0){
		tribipyrStructHg::gen_summary_report_tribipyr_hg(no_of_compared_tribipyr_hg);
	}

	if (no_of_compared_tribipyr_cd!=0){
		tribipyrStructCd::gen_summary_report_tribipyr_cd(no_of_compared_tribipyr_cd);
	}

	if (no_of_compared_tribipyr_ni!=0){
		tribipyrStructNi::gen_summary_report_tribipyr_ni(no_of_compared_tribipyr_ni);
	}
}



