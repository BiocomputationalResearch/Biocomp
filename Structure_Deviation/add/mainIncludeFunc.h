
#ifndef MAININCLUDEFUNC_H_
#define MAININCLUDEFUNC_H_
#include "coordinationsite.h"

bool testNumber(char *str);
string convert_to_mol2(string, string, string);
void form_tripl_sites(vector<metal_binding_site>, string, string);
void form_tetrapl_sites(vector<metal_binding_site>, string, string);
void form_tetra_sites(vector<metal_binding_site>, string, string);
void form_tribipyr_sites(vector<metal_binding_site>, string, string);
void form_tetrapyr_sites(vector<metal_binding_site>, string, string);
void form_octa_sites(vector<metal_binding_site>, string, string);
void generate_summary_report_tripl();
void generate_summary_report_tetrapl();
void generate_summary_report_tetra();
void generate_summary_report_tribipyr();
void generate_summary_report_tetrapyr();
void generate_summary_report_octa();


#endif /* MAININCLUDEFUNC_H_ */
