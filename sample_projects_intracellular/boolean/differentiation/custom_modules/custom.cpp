/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/
#define _USE_MATH_DEFINES
#include <cmath>
#include <sstream>
#include "./custom.h"
#include "../BioFVM/BioFVM.h"  

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  

	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	

	return; 
}

void setup_tissue( void )
{
	// place a cluster of tumor cells at the center 
 	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	double tumor_radius = parameters.doubles("tumor_radius");
	// Parameter<double> temp; 
	
	Cell* pCell = NULL; 
	
	std::vector<std::vector<double>> positions;	

	if (default_microenvironment_options.simulate_2D == true){
		positions = create_cell_disc_positions(cell_radius,tumor_radius);
		std::cout << "ENABLED 2D SIMULATION" << std::endl; 
	}
	else
		positions = create_cell_sphere_positions(cell_radius,tumor_radius);
	std::cout << "creating " << positions.size() << " closely-packed tumor cells ... " << std::endl;


	for( int i=0; i < positions.size(); i++ )
	{
		pCell = create_cell(get_cell_definition("cancercell")); 
		
		pCell->assign_position( positions[i] );
	} 

	load_cells_from_pugixml();
	
	return; 
}


std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius)
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*2;
	double z_spacing= cell_radius*sqrt(3);

	std::vector<double> tempPoint(3,0.0);
	// std::vector<double> cylinder_center(3,0.0);

	for(double z=-sphere_radius;z<sphere_radius;z+=z_spacing, zc++)
	{
		for(double x=-sphere_radius;x<sphere_radius;x+=x_spacing, xc++)
		{
			for(double y=-sphere_radius;y<sphere_radius;y+=y_spacing, yc++)
			{
				tempPoint[0]=x + (zc%2) * 0.5 * cell_radius;
				tempPoint[1]=y + (xc%2) * cell_radius;
				tempPoint[2]=z;

				if(sqrt(norm_squared(tempPoint))< sphere_radius)
				{ cells.push_back(tempPoint); }
			}

		}
	}
	return cells;

}

std::vector<std::vector<double>> create_cell_disc_positions(double cell_radius, double disc_radius)
{	 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	
	double x = 0.0; 
	double y = 0.0; 
	double x_outer = 0.0;

	std::vector<std::vector<double>> positions;
	std::vector<double> tempPoint(3,0.0);
	
	int n = 0; 
	while( y < disc_radius )
	{
		x = 0.0; 
		if( n % 2 == 1 )
		{ x = 0.5 * cell_spacing; }
		x_outer = sqrt( disc_radius*disc_radius - y*y ); 
		
		while( x < x_outer )
		{
			tempPoint[0]= x; tempPoint[1]= y;	tempPoint[2]= 0.0;
			positions.push_back(tempPoint);			
			if( fabs( y ) > 0.01 )
			{
				tempPoint[0]= x; tempPoint[1]= -y;	tempPoint[2]= 0.0;
				positions.push_back(tempPoint);
			}
			if( fabs( x ) > 0.01 )
			{ 
				tempPoint[0]= -x; tempPoint[1]= y;	tempPoint[2]= 0.0;
				positions.push_back(tempPoint);
				if( fabs( y ) > 0.01 )
				{
					tempPoint[0]= -x; tempPoint[1]= -y;	tempPoint[2]= 0.0;
					positions.push_back(tempPoint);
				}
			}
			x += cell_spacing; 
		}		
		y += cell_spacing * sqrt(3.0)/2.0; 
		n++; 
	}
	return positions;
}

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ 	
	return; 	
} 

void pre_update_intracellular(Cell* pCell, Phenotype& phenotype, double dt){

	return;
}

void post_update_intracellular(Cell* pCell, Phenotype& phenotype, double dt){


	return;
}

	// FUNCTIONS TO PLOT CELLS

std::vector<std::string> ECM_coloring_function( Cell* pCell)
{
	std::vector< std::string > output( 4 , "black" );
	std::string parameter = parameters.strings("parameter_to_visualize");;
	double param = pCell->custom_data[parameter];
	int color = (int) round( param * 255.0 );
	if(color > 255){
		color = 255;
	}
	char szTempString [128];
	sprintf( szTempString , "rgb(%u,0,%u)", color, 255-color );
	output[0].assign( szTempString );
	return output;

}

std::vector<std::string> phase_coloring_function( Cell* pCell )
{
	std::vector< std::string > output( 4 , "rgb(0,0,0)" );

	if ( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::Ki67_negative )
	{
		output[0] = "rgb(0,255,0)"; //green
		output[2] = "rgb(0,125,0)";
		
	}
	if ( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::Ki67_positive_premitotic )
	{
		output[0] = "rgb(255,0,0)"; //red
		output[2] = "rgb(125,0,0)";
		
	}
	if ( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::Ki67_positive_postmitotic )
	{
		output[0] = "rgb(255,255,0)"; //yellow
		output[2] = "rgb(125,125,0)";
		
	}
	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
	{
		output[0] = "rgb(165,42,42)"; //brown
		output[2] = "rgb(165,42,42)";
	}
	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic )  
	{
		output[0] = "rgb(0,0,0)";
		output[2] = "rgb(0,0,0)";
	}
	return output;
}

std::vector<std::string> node_coloring_function( Cell* pCell )
{
	std::vector< std::string > output( 4 , "rgb(0,0,0)" );
	//std::cout << pCell->phenotype.intracellular->get_boolean_node_value( parameters.strings("node_to_visualize"));
	if ( !pCell->phenotype.intracellular->get_boolean_variable_value( parameters.strings("node_to_visualize") ) ) //node off
	{
		output[0] = "rgb(0,0,255)"; //blue
		output[2] = "rgb(0,0,125)";
		
	}
	else{
		output[0] = "rgb(0,255,0)"; //green
		output[2] = "rgb(0,125,0)";
	}
	
	return output;
}


std::vector<std::string> my_coloring_function( Cell* pCell )
{
	
	 int color_number = parameters.ints("color_function");

	 if (color_number == 0)
	 	return ECM_coloring_function(pCell);
	 if (color_number == 1)
	 	return phase_coloring_function(pCell);
	 else 
	 	return node_coloring_function( pCell );
}

std::vector<std::string> my_coloring_function_for_stroma( double concentration, double max_conc, double min_conc )
{
	 return paint_by_density_percentage( concentration,  max_conc,  min_conc); 

}
