// ======================================================================
//
// file          : pbes2bes
// date          : 24-02-2007
// version       : 0.1.0
//
// author(s)     : Alexander van Dam <avandam@damdonk.nl>
//
// ======================================================================


#define NAME "pbes2bes"
#define VERSION "0.1.0"
#define AUTHOR "Alexander van Dam"


//C++
#include <cstdio>
#include <exception>
#include <iostream>
#include <fstream>
#include <string>
#include <utility>

#include <sstream>

//Boost
#include <boost/program_options.hpp>

//MCRL-specific
#include "liblowlevel.h"
#include "libprint_c.h"

//LPS-Framework
#include "lps/pbes.h"
#include "lps/pbes_utility.h"
#include "lps/data_operators.h"
#include "lps/sort.h"
//#include "lps/sort_utility.h"

//ATERM-specific
#include "atermpp/substitute.h"
#include "lps/identifier_string.h"
#include "atermpp/utility.h"
#include "atermpp/indexed_set.h"
#include "atermpp/table.h"

//Tool-specific
#include "pbes_rewrite.h"		// PBES rewriter
#include "sort_functions.h"

using namespace std;
using namespace lps;

using atermpp::make_substitution;

namespace po = boost::program_options;

//Type definitions
//----------------
typedef struct{
	string opt_outputformat;
	string infilename;
	string outfilename;
} t_tool_options;

typedef struct {
	data_variable_list finite_var;		// List of all finite variables
	data_variable_list infinite_var;	// List of all infinite variables
	data_expression_list finite_exp;	// List of all finite expressions
	data_expression_list infinite_exp;	// List of all infinite expressions
} t_instantiations;

//Function declarations used by main program
//------------------------------------------
static t_tool_options parse_command_line(int argc, char** argv);
//Post: The command line options are parsed.
//      The program has aborted with a suitable error code, if:
//		- errors were encounterd
//		- non-standard behaviour was requested (help, version)
//Ret:	The parsed command line options

static pbes create_bes(pbes pbes_spec, t_tool_options tool_options);
//Post: tool_options.infilename contains a PBES ("-" indicates stdin)
//Ret:	The BES generated from the PBES

pbes load_pbes(t_tool_options tool_options);
//Post: tool_options.infilename contains a PBES ("-" indicates stdin)
//Ret: The pbes loaded from infile

void save_pbes(pbes pbes_spec, t_tool_options tool_options);
//Post: pbes_spec is saved
//Ret: -

void save_pbes_in_cwi_format(pbes pbes_spec, string outfilename);
//Post: pbes_spec is saved in cwi-format
//Ret: -

string convert_rhs_to_cwi(pbes_expression p, atermpp::indexed_set *variables);
// Function used to convert a pbes_expression to the variant used by the cwi-output

pbes do_naive_algorithm(pbes pbes_spec, t_tool_options tool_options);
//

propositional_variable_instantiation create_naive_propositional_variable_instantiation(propositional_variable_instantiation propvarinst, atermpp::table *enumerated_sorts);
// Create a propositional variable instantiation with the checks needed in the naive algorithm

identifier_string create_propvar_name(identifier_string propvar_name, data_expression_list finite_exp);
// Create a identifier string containing the name for use in the propositional_variable(_instantiation)

//Main Program
//------------
int main(int argc, char** argv)
{
	//Initialise ATerm library and lowlevel-functions
	ATerm bottom;
	ATinit(argc, argv, &bottom);
	gsEnableConstructorFunctions();

	//Parse command line
	t_tool_options tool_options = parse_command_line(argc, argv);

	//Load PBES
	pbes pbes_spec = load_pbes(tool_options);

	//Process the pbes
	pbes_spec = create_bes(pbes_spec, tool_options);

	//Save PBES
	save_pbes(pbes_spec, tool_options);

	return 0;
}

//function create_bes
//-------------------
pbes create_bes(pbes pbes_spec, t_tool_options tool_options)
{
	//Do naive algorithm. In later stage a check for naive or improved is added.
	if (!pbes_spec.is_well_formed())
	{
		gsErrorMsg("The PBES is not well formed. Pbes2bes cannot handle this kind of PBES's\nComputation aborted.");
		exit(1);
	}
	if (!pbes_spec.is_closed())
	{
		gsErrorMsg("The PBES is not closed. Pbes2bes cannot handle this kind of PBES's\nComputation aborted.");
		exit(1);
	}
	pbes_spec = do_naive_algorithm(pbes_spec, tool_options);

	//return new pbes
	return pbes_spec;
}

//function do_naive_algorithm
//---------------------------
pbes do_naive_algorithm(pbes pbes_spec, t_tool_options tool_options)
{
	// Get all parts from the PBES (data specification, equation system, initial state)
	propositional_variable_instantiation initial_state = pbes_spec.initial_state();
	equation_system eqsys = pbes_spec.equations();
	data_specification data = pbes_spec.data();

	// Needed variables in the whole function
	equation_system result_eqsys;				// resulting equation system
	bool is_bes_system = true;					// True if resulting eqsys is a BES
	int nr_of_equations = 0;					// Nr of equations computed
	Rewriter *rewriter = createRewriter(data); 	// Data rewriter

	// Empty data_variable_list and data_expression_list
	data_variable_list empty_data_variable_list;
	data_expression_list empty_data_expression_list;

	atermpp::table *sort_enumerations = new atermpp::table(10,50);

	//Populate sort_enumerations with all enumerations for finite sorts
	for (equation_system::iterator eq_i = eqsys.begin(); eq_i != eqsys.end(); eq_i++)
	{
		data_variable_list parameters = eq_i->variable().parameters();
		for (data_variable_list::iterator p = parameters.begin(); p != parameters.end(); p++)
		{
			lps::sort current_sort = p->sort();
			if (sort_enumerations->get(current_sort) == NULL)
			{ // If the sort is not in the enumerated_sorts
				if (check_finite(data.constructors(), current_sort))
				{ // If the sort is finite
					data_expression_list enumerations_from_sort = enumerate_constructors(data.constructors(), current_sort);
					sort_enumerations->put(current_sort, enumerations_from_sort);
				}
			}
		}
	}

	// Create the rewritten equation system
	for (equation_system::iterator eq_i = eqsys.begin(); eq_i != eqsys.end(); eq_i++)
	{
		// Get current equation
		pbes_equation equation = *eq_i;
		// Get all parts from the PBES_EQUATION (propositional_variable, formula)
		propositional_variable propvar = equation.variable();
		pbes_expression formula = equation.formula();

		// Get the name and parameters of the propvar
		identifier_string propvar_name = propvar.name();
		data_variable_list propvar_parameters = propvar.parameters();

		// Verbose print the equation which is currently done
		string propvar_name_string = propvar_name;
		gsVerboseMsg("Rewriting PBES equation with propositional variable %s\n", propvar_name_string.c_str());

		// Vector of instantiations
		vector< t_instantiations > instantiation_list;
		t_instantiations instantiation;			// Needed?
		t_instantiations current_values;		// Current results

		// Add empty instantiation to the list
		instantiation_list.push_back(instantiation);

		// For each variable of the equation
		for (data_variable_list::iterator p = propvar_parameters.begin(); p != propvar_parameters.end(); p++)
		{
			// Vector of instantiations for intermediate results
			vector< t_instantiations > intermediate_instantiation_list;
			if (sort_enumerations->get(p->sort()) == NULL)
			{ // The sort is infinite
				is_bes_system = false;		// Resulting system is not a BES
				current_values.infinite_var = push_back(current_values.infinite_var, *p);
				for (vector< t_instantiations >::iterator inst_i = instantiation_list.begin(); inst_i != instantiation_list.end(); inst_i++)
				{
					// Add current_values to intermediate_instantiation_list
					intermediate_instantiation_list.push_back(current_values);
				}
			}
			else
			{ // The sort is finite
				current_values.finite_var = push_back(current_values.finite_var, *p);
				data_expression_list enumerations = sort_enumerations->get(p->sort());

				for (vector< t_instantiations >::iterator inst_i = instantiation_list.begin(); inst_i != instantiation_list.end(); inst_i++)
				{
					for (data_expression_list::iterator e = enumerations.begin(); e != enumerations.end(); e++)
					{
						current_values.finite_exp = push_back(inst_i->finite_exp, *e);
						intermediate_instantiation_list.push_back(current_values);
					}
				}
			}
			instantiation_list = intermediate_instantiation_list;
		}

		// All instantiations for the current equation are computed. Now make equations out of them
		for (vector< t_instantiations >::iterator inst_i = instantiation_list.begin(); inst_i != instantiation_list.end(); inst_i++)
		{
			// Create the needed propvar
			propositional_variable propvar_current = propositional_variable(create_propvar_name(propvar_name, inst_i->finite_exp), inst_i->infinite_var);

			pbes_expression current_expression = formula;	// Current expression
			// Substitute all variables which are instantiated
			current_expression = current_expression.substitute(make_list_substitution(inst_i->finite_var, inst_i->finite_exp));
			// Rewrite the rhs as far as possible
			current_expression = pbes_expression_rewrite(current_expression, data, rewriter);
			// Lists to replace variables in the rhs
			propositional_variable_instantiation_list oldpropvarinst_list;
			propositional_variable_instantiation_list newpropvarinst_list;
			// Get all propvarinst of the rhs
			set< propositional_variable_instantiation > propvarinst_set = find_propositional_variable_instantiations(current_expression);

			// For each propvarinst in the set
			for (set< propositional_variable_instantiation >::iterator pvi = propvarinst_set.begin(); pvi != propvarinst_set.end(); pvi++)
			{
				// Put the propvarinst in the old list
				oldpropvarinst_list = push_front(oldpropvarinst_list, *pvi);
				// Create the new propvarinst
				propositional_variable_instantiation newpropvarinst = create_naive_propositional_variable_instantiation(*pvi, sort_enumerations);
				// Add the new propvarinst in the new list
				newpropvarinst_list = push_front(newpropvarinst_list, newpropvarinst);
			}

			// Replace the propvarinsts with the new ones
			current_expression = current_expression.substitute(make_list_substitution(oldpropvarinst_list, newpropvarinst_list));

			// Add the new pequation tin the resulting system
			result_eqsys.push_back(pbes_equation(eq_i->symbol(), propvar_current, current_expression));

			// Verbose messages after each 100 equations
			if (++nr_of_equations % 100 == 0)
				gsVerboseMsg("At Boolean equation %d\n", nr_of_equations);

		}
	}

	// rewrite initial state
	propositional_variable_instantiation new_initial_state = create_naive_propositional_variable_instantiation(initial_state, sort_enumerations);

	// create resulting pbes
	pbes result = pbes(data, result_eqsys, new_initial_state);

	gsVerboseMsg("Number of Boolean equation systems computed: %d\n", nr_of_equations);

	delete rewriter;

	return result;
}


//function create_naive_propositional_variable_instantiation
//----------------------------------------------------------
propositional_variable_instantiation create_naive_propositional_variable_instantiation(propositional_variable_instantiation propvarinst, atermpp::table *enumerated_sorts)
{
	data_expression_list finite_expression;
	data_expression_list infinite_expression;
	//For each parameter:
	for (data_expression_list::iterator p = propvarinst.parameters().begin(); p != propvarinst.parameters().end(); p++)
	{
		if (enumerated_sorts->get(p->sort()) != NULL)
		{
			//If it is of finite sort, add it to the list of finite expressions of the propvarinst
			if (is_data_operation(*p))
			{ // If p is a correct data operation
				finite_expression = push_back(finite_expression, *p);
			}
			else if (is_data_variable(*p))
			{ // If p is a freevar
				gsErrorMsg("The PBES contains a free variable in the propositional variable instantiation.\n");
				gsErrorMsg("Try using mcrl22lps with the flag -f or --no-freevars to solve this.\n");
				exit(1);
			}
		}
		else
		{
			//If it is of infinite sort, add it to the list of infinite expressions of the propvarinst
			infinite_expression = push_back(infinite_expression, *p);
		}
	}

	return propositional_variable_instantiation(create_propvar_name(propvarinst.name(), finite_expression), infinite_expression);
}

//function create_propvar_name
//----------------------------
identifier_string create_propvar_name(identifier_string propvar_name, data_expression_list finite_exp)
{
	// string representation of the original propvar_name
	string propvar_name_current = propvar_name;
	// Only add data to the string if the finite_exp list is non-empty. Currently "_" is used as seperator
	if (!finite_exp.empty())
	{
		// Add each parameter
		for (data_expression_list::iterator del_i = finite_exp.begin(); del_i != finite_exp.end(); del_i++)
		{
			propvar_name_current += "_";
			propvar_name_current += data_operation(del_i->head()).name();
		}
	}

	return propvar_name_current;
}

//function load_pbes
//------------------
pbes load_pbes(t_tool_options tool_options)
{
	string infilename = tool_options.infilename;

	pbes pbes_spec;
	if (infilename == "-")
	{
		if (!pbes_spec.load("-"))
		{
			gsErrorMsg("Cannot open PBES from stdin\n");
			exit(1);
		}
	}
	else
	{
		if (!pbes_spec.load(infilename))
		{
			gsErrorMsg("Cannot open PBES from '%s'\n", infilename.c_str());
			exit(1);
		}
	}
	return pbes_spec;
}

//function save_pbes
//------------------
void save_pbes(pbes pbes_spec, t_tool_options tool_options)
{
	//Outfile string
	string outfilename = tool_options.outfilename;

	//Write PBES
	if (tool_options.opt_outputformat == "external")
	{ // In external format
		if (!pbes_spec.save(outfilename, false))
		{
			gsErrorMsg("Could not save PBES to %s\n", outfilename.c_str());
			exit(1);
		}
	}
	else if (tool_options.opt_outputformat == "binary")
	{ // in binary format
		if (!pbes_spec.save(outfilename, true))
		{
			gsErrorMsg("Could not save PBES to %s\n", outfilename.c_str());
			exit(1);
		}
	}
	else if (tool_options.opt_outputformat == "cwi")
	{ // in CWI format only if the result is a BES, otherwise Binary
		if (pbes_spec.equations().is_bes())
		{
			save_pbes_in_cwi_format(pbes_spec, outfilename);
		}
		else
		{
			gsWarningMsg("Result is not a BES. Result is written in binary format.\n");
			if (!pbes_spec.save(outfilename, true))
			{
				gsErrorMsg("Could not save PBES to %s\n", outfilename.c_str());
				exit(1);
			}

		}
	}
}

//function save_pbes_in_cwi_format
//--------------------------------
void save_pbes_in_cwi_format(pbes pbes_spec, string outfilename)
{
	// Use an indexed set to keep track of the variables and their cwi-representations
	equation_system eqsys = pbes_spec.equations();
	atermpp::indexed_set *variables = new atermpp::indexed_set(2*eqsys.size(), 50);
	// Fill the indexed set
	for (equation_system::iterator eq = eqsys.begin(); eq != eqsys.end(); eq++)
	{
		variables->put(eq->variable());
	} // because pbes is closed, all variables are done at this point

	ofstream outputfile;
   	outputfile.open(outfilename.c_str(), ios::trunc);
	// Rewrite all equations to CWI format
	for (equation_system::iterator eq = eqsys.begin(); eq != eqsys.end(); eq++)
	{
		// fixpoint:
		// 	mu => min
		// 	nu => max
		string fp;
		(eq->symbol().is_mu())?fp = "min":fp = "max";
		// Check if the
		//if (eq->variable().parameters().size() != 0)
		//{
		//	gsErrorMsg("The used equation system is not a BES. Could not save this in CWI-format.\n");
		//	exit(1);
		//}

		// variable:
		// 	Integer representation of the propositional variable of the left hand side
	    string variable;
		stringstream variable_stream;
		variable_stream << variables->index(eq->variable());
		variable = variable_stream.str();

		// Convert the right hand side of the equation to the CWI format
		string rhs = convert_rhs_to_cwi(eq->formula(), variables);

		// Create the equation
		string equation = fp + " " + variable + " = " + rhs + "\n";

		if (outputfile.is_open())
			outputfile << equation;
		else
		{
			gsErrorMsg("Could not save PBES to %s\n", outfilename.c_str());
			exit(1);
		}
	}
	outputfile.close();
}

//function convert_rhs_to_cwi
//---------------------------
string convert_rhs_to_cwi(pbes_expression p, atermpp::indexed_set *variables)
{
	string result;
	if (pbes_expr::is_true(p))
		// PBESTrue => T
		result = "T";
	else if (pbes_expr::is_false(p))
		// PBESFalse => F
		result = "F";
	else if (pbes_expr::is_and(p))
	{
		//PBESAnd(a,b) => (a & b)
		string left = convert_rhs_to_cwi(pbes_expr::lhs(p), variables);
		string right = convert_rhs_to_cwi(pbes_expr::rhs(p), variables);
/*		if ((left == "false") || (right == "false"))
			result = "false";
		else if (left == "true")
			result = right;
		else if (right == "true")
			result = left;
		else
*/
		result = "(" + left + "&" + right + ")";
	}
	else if (pbes_expr::is_or(p))
	{
		//PBESOr(a,b) => (a | b)
		string left = convert_rhs_to_cwi(pbes_expr::lhs(p), variables);
		string right = convert_rhs_to_cwi(pbes_expr::rhs(p), variables);
/*		if ((left == "true") || (right == "true"))
			result = "true";
		else if (left == "false")
			result = right;
		else if (right == "false")
			result = left;
		else
*/
		result = "(" + left + "|" + right + ")";
	}
	else if (pbes_expr::is_propositional_variable_instantiation(p))
	{
		// PropVar => <Int>
		propositional_variable_instantiation propvarinst = propositional_variable_instantiation(p);
		data_variable_list empty;
		propositional_variable propvar = propositional_variable(propvarinst.name(), empty);
		long variable = variables->index(propvar);
		if (variable < 0)
		{
			gsErrorMsg("Error: The BES is not closed. Write to cwi-format failed.");
			exit(1);
		}
		else
		{
			stringstream variable_stream;
			variable_stream << variable;
			result = variable_stream.str();
		}
	}
	else
	{
		gsErrorMsg("The used equation system is not a BES. Could not save this in CWI-format.\n");
		exit(1);
	}
	return result;
}

//function parse_command_line
//---------------------------
t_tool_options parse_command_line(int argc, char** argv)
{
	t_tool_options tool_options;
	string opt_outputformat;
	vector< string > file_names;

	po::options_description desc;

	desc.add_options()
			("output,o",	po::value<string>(&opt_outputformat)->default_value("binary"), "use outputformat arg (default 'binary');"
			 				"available outputformats are binary, external and cwi")
			("verbose,v",	"turn on the display of short intermediate messages")
			("debug,d",		"turn on the display of detailed intermediate messages")
			("version",		"display version information")
			("help,h",		"display this help")
			;

	po::options_description hidden("Hidden options");
	hidden.add_options()
			("file_names",		po::value< vector< string > >(), "input/output files")
			;

	po::options_description cmdline_options;
	cmdline_options.add(desc).add(hidden);

	po::options_description visible("Allowed options");
	visible.add(desc);

	po::positional_options_description p;
	p.add("file_names", -1);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("help"))
	{
		cerr << "Usage: " << argv[0] << " [OPTION]... [INFILE [OUTFILE]]" << endl;
		cerr << "Rewrites PBES from stdin or INFILE to an equivalent BES." << endl;
		cerr << "By default the result is written to stdout, otherwise to OUTFILE." << endl;
		cerr << endl;
		cerr << desc;

		exit(0);
	}

	if (vm.count("version"))
	{
		cerr << NAME << " " << VERSION <<  " (revision " << REVISION << ")" << endl;

		exit(0);
	}

	if (vm.count("debug"))
	{
		gsSetDebugMsg();
	}

	if (vm.count("verbose"))
	{
		gsSetVerboseMsg();
	}

	if (vm.count("output")) // Output format
	{
		opt_outputformat = vm["output"].as< string >();
		if (!((opt_outputformat == "external") || (opt_outputformat == "binary") || (opt_outputformat == "cwi")))
		{
			gsErrorMsg("Unknown outputformat specified. Available outputformats are binary, external and cwi");
			exit(1);
		}
	}

	if (vm.count("file_names"))
	{
		file_names = vm["file_names"].as< vector< string > >();
	}

	// Check for number of arguments and infilename/outfile
	string infilename;
	string outfilename;
	if (file_names.size() == 0)
	{
		/* Read from standard input */
		infilename = "-";
		outfilename = "-";
	}
	else if (2 < file_names.size())
	{
		cerr << NAME <<": Too many arguments" << endl;
		exit(1);
	}
	else
	{
		infilename = file_names[0];
		if (file_names.size() == 2)
		{
			outfilename = file_names[1];
		}
		else
		{
			outfilename = "-";
		}
	}
	tool_options.opt_outputformat = opt_outputformat;
	tool_options.infilename = infilename;
	tool_options.outfilename = outfilename;
	return tool_options;
}


