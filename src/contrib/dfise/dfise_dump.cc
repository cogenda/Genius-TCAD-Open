/********************************************************************************/
/*     888888    888888888   88     888  88888   888      888    88888888       */
/*   8       8   8           8 8     8     8      8        8    8               */
/*  8            8           8  8    8     8      8        8    8               */
/*  8            888888888   8   8   8     8      8        8     8888888        */
/*  8      8888  8           8    8  8     8      8        8            8       */
/*   8       8   8           8     8 8     8      8        8            8       */
/*     888888    888888888  888     88   88888     88888888     88888888        */
/*                                                                              */
/*       A Three-Dimensional General Purpose Semiconductor Simulator.           */
/*                                                                              */
/*                                                                              */
/*  Copyright (C) 2007-2008                                                     */
/*  Cogenda Pte Ltd                                                             */
/*                                                                              */
/*  Please contact Cogenda Pte Ltd for license information                      */
/*                                                                              */
/*  Author: Gong Ding   gdiso@ustc.edu                                          */
/*                                                                              */
/********************************************************************************/
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <fstream>


#include "vector_value.h"
#include "tensor_value.h"
#include "dfise.h"

#ifdef WINDOWS
  #include <io.h>      // for windows _access function
#else
  #include <unistd.h>  // for POSIX access function
#endif


void printusage()
{
  std::cout<<"Usage: dump_dfise [-i input] [-l] [-o output] [-v variable]\n";
  std::cout<<"Options\n";
  std::cout<<"  -h\t\tDisplay this help\n";
  std::cout<<"  -i\t\tInput DFISE file without extersion (xxx for xxx.grd and xxx.dat)\n";
  std::cout<<"  -o\t\tOutput data file, default is result.dat\n";
  std::cout<<"  -l\t\tList all the variables in DFISE file\n";
  std::cout<<"  -v\t\tSpecify variable index for extracting\n";
  std::cout<<"Contact support@cogenda.com for more helps.\n";
}

/**
 * a small tool for dump data in dfise file
 */
int main(int argc, char **argv)
{

  std::string input_file;
  std::string output_file("result.dat");

  int variable_index = -1;
  bool list_variables=false;

  std::cout<<"*********************************************************************\n";
  std::cout<<"*  dump_dfise is a tool for extracting data in Synopsys DFISE file. *\n";
  std::cout<<"*  Copyright (C) 2009-2010 by Cogenda EDA.                          *\n";
  std::cout<<"*              http://www.cogenda.com/                              *\n";
  std::cout<<"*********************************************************************\n";

  if(argc == 1)
  {
    printusage();
    exit(1);
  }

  // parse command line
  for(int i=1; i<argc; i++)
  {
    if(argv[i][0] == '-')
    {
      // input file
      if(strncmp(argv[i], "-i", strlen("-i"))==0)
      {
        if(i==argc-1)
        {
          std::cerr<<"Error: -i switch given but no input file specified.\n\n";
          printusage();
          exit(1);
        }
        else
        {
          input_file = argv[i+1];
        }
      }

      // output file
      if(strncmp(argv[i], "-o", strlen("-o"))==0)
      {
        if(i==argc-1)
        {
          std::cerr<<"Error: -o switch given but no output file specified.\n\n";
          printusage();
          exit(1);
        }
        else
        {
          output_file = argv[i+1];
        }
      }

      if(strncmp(argv[i], "-l", strlen("-l"))==0)
      {
        list_variables = true;
      }

      if(strncmp(argv[i], "-v", strlen("-v"))==0)
      {
        if(i==argc-1)
        {
          std::cerr<<"Error: -v switch given but no variable index specified.\n\n";
          printusage();
          exit(1);
        }
        else
        {
          std::stringstream ss;
          ss << argv[i+1];
          ss >> variable_index;
        }
      }

      // help info
      if(strncmp(argv[i], "-h", strlen("-h"))==0)
      {
        printusage();
        exit(1);
      }

    }
  } // end of parse command line

  {
    std::string grid_file = input_file + ".grd";
    std::string data_file = input_file + ".dat";
#ifdef WINDOWS
    if ( _access( (void*)grid_file.c_str(),  04 ) == -1 )
#else
    if ( access( grid_file.c_str(),  R_OK ) == -1 )
#endif
    {
      std::cerr<<"ERROR: I can't read DFISE Grid file "<<grid_file <<", access failed.\n";
      exit(0);
    }

#ifdef WINDOWS
    if ( _access( (void*)data_file.c_str(),  04 ) == -1 )
#else
    if ( access( data_file.c_str(),  R_OK ) == -1 )
#endif
    {
      std::cerr<<"ERROR: I can't read DFISE Data file "<<data_file <<", access failed.\n";
      exit(0);
    }
  }

  // parse dfise file
  std::cout<<"Parsing DFISE file:\n";
  DFISE::DFISE_MESH ise_reader;
  ise_reader.parse_dfise(input_file);

  const DFISE::INFO & grid_info = ise_reader.get_grid_info();
  const DFISE::GRID & grid      = ise_reader.get_grid();

  bool no_variable_given =  variable_index < 0;

  // output all the variables
  if(list_variables || no_variable_given)
  {
    std::cout<<"List all the variables in this DFISE file:\n";
    unsigned int n_datasets = ise_reader.n_datasets();

    //for each dataset
    for(unsigned int n=0; n<n_datasets; ++n)
    {
      const DFISE::DATASET * dataset = ise_reader.get_dataset(n);
      std::cout << std::setw(3) << n << " " << dataset->name << std::endl;
      std::cout << "      Defined in region(s):" << std::endl;
      for(unsigned int r=0; r<dataset->validity.size(); ++r)
      {
        std::cout << "        " << dataset->validity[r] << std::endl;
      }
      std::cout << std::endl;
    }
  }

  // input variable by user
  if( no_variable_given )
  {
    std::cout<<"Which variable to dump [0-"<<ise_reader.n_datasets()-1<<"]:";
    std::cin>>variable_index;
  }


  if( variable_index<0 || variable_index >= ise_reader.n_datasets() )
  {
    std::cerr<<"ERROR: Variable index "<< variable_index <<" out of range.\n";
    exit(0);
  }

  // write variable data at its location
  const DFISE::DATASET * dataset = ise_reader.get_dataset(variable_index);
  std::string variable = dataset->name;

  std::cout<<"Writing " << variable << " to file " << output_file <<"..."<< std::endl;

  std::ofstream fout;
  fout.open ( output_file.c_str(), std::ofstream::trunc );
  fout<<std::scientific;

  // this map hold the node index to value map, this is what we wnat
  VectorValue<double> translate(grid.translate);
  TensorValue<double> transform(grid.transform);
  const std::map<unsigned int, unsigned int> & node_to_value_index_map = dataset->node_to_value_index_map;
  std::map<unsigned int, unsigned int>::const_iterator it = node_to_value_index_map.begin();
  for(; it!=node_to_value_index_map.end(); ++it)
  {
    const DFISE::Point & point = grid.Vertices[it->first];
    VectorValue<double> p(point[0], point[1], point[2]);
    p = transform*p + translate;
    double value = dataset->get_scaler_value_by_node_index(it->first);
    fout << p(0) << std::setw(15)
         << p(1) << std::setw(15)
         << p(2) << std::setw(15)
         << value
         << std::endl;
  }

  fout.close();



}
