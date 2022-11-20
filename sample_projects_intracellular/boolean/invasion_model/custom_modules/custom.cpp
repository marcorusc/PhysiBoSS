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

	cell_defaults.functions.pre_update_intracellular = pre_update_intracellular;
	cell_defaults.functions.post_update_intracellular = post_update_intracellular;
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 

	Cell_Definition* pCD = find_cell_definition( "epithelial");

	pCD->functions.pre_update_intracellular = pre_update_intracellular;
	pCD->functions.post_update_intracellular = post_update_intracellular;
	pCD->functions.custom_cell_rule = custom_function; 
	pCD->functions.contact_function = contact_function; 
	pCD->functions.update_velocity = custom_update_cell_velocity; 

	pCD = find_cell_definition( "mesenchymal");
	pCD->functions.pre_update_intracellular = pre_update_intracellular;
	pCD->functions.post_update_intracellular = post_update_intracellular;
	pCD->functions.custom_cell_rule = custom_function; 
	pCD->functions.contact_function = contact_function;
	pCD->functions.update_velocity = custom_update_cell_velocity; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void set_substrate_density(int density_index, double max, double min, double radius)
{
	std::cout << "SETTING SUBSTRATE --> " << density_index << std::endl;
	// Inject given concentration on the extremities only

	std::cout << microenvironment.number_of_voxels() << "\n";

	for (unsigned int n = 0; n < microenvironment.number_of_voxels(); n++)
	{
		auto current_voxel = microenvironment.voxels(n);
		double t_norm = norm(current_voxel.center);

		if ((radius - t_norm) <= 0)
			microenvironment.density_vector(n)[density_index] = current_value(min, max, uniform_random());
	}
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 

	double radius_ECM = parameters.doubles("config_radius");
	double radius_tgfbeta = parameters.doubles("tgfbeta_radius");

	double ECM_min = parameters.doubles("density_ECM_min");
	double ECM_max = parameters.doubles("density_ECM_max");
	double tgfbeta_max = parameters.doubles("density_tgfbeta_max");
	double tgfbeta_min = parameters.doubles("density_tgfbeta_min");

	int ecm_index = microenvironment.find_density_index("ecm");
	int tgfbeta_index = microenvironment.find_density_index("tgfbeta");
	
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	set_substrate_density(ecm_index, ECM_max, ECM_min, radius_ECM);
	
	//set_substrate_density(tgfbeta_index, tgfbeta_max, tgfbeta_min, radius_tgfbeta);

	return; 
}

void setup_tissue( void )
{
	// place a cluster of tumor cells at the center 
	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	
	double tumor_radius = parameters.doubles( "tumor_radius" ); // 250.0; 
	
	// Parameter<double> temp; 
	
	int i = parameters.doubles.find_index( "tumor_radius" ); 
	
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
		pCell = create_cell(get_cell_definition("epithelial")); 
		
		pCell->assign_position( positions[i] );
	}
	
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


//std::vector<std::string> my_coloring_function( Cell* pCell )
//{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ 	
	return; 	
} 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ 
	std::vector<double> displacement = pOther->position;
	displacement -= pMe->position;
	double distance = norm( displacement ); 
			
	double max_distance = pMe->phenotype.geometry.radius + 
				pOther->phenotype.geometry.radius; 
	max_distance *=  1.1;  //parameters.doubles("max_interaction_factor"); 

			//std::cout << max_distance << " - " << distance << "\n";

	double interaction_distance = max_distance - distance;
	if (interaction_distance < 0){
		double perc_distance = distance / pMe->phenotype.geometry.radius ;
		pMe->custom_data["cell_contact"] += perc_distance;
			}

	return; 
} 

void pre_update_intracellular(Cell* pCell, Phenotype& phenotype, double dt){

	//std::cout << "how much ? " << pCell->custom_data["ecm_contact"] << std::endl;

	return;
}

void post_update_intracellular(Cell* pCell, Phenotype& phenotype, double dt){

	// to insert trigger of the degradation according to the value of ecm_degrad

	if ( pCell->phenotype.intracellular->has_variable("ECM_degrad") ){
	 	set_mmp(pCell, pCell->phenotype.intracellular->get_boolean_variable_value("ECM_degrad") );}

	// to insert trigger of the migration according to the value of the node Migration

	if ( pCell->phenotype.intracellular->has_variable("Migration"))
	{
		set_oxygen_motility(pCell, pCell->phenotype.intracellular->get_boolean_variable_value("Migration"));
	}

	// to insert trigger of the detatch/attach of the cell according to the value of node EMT OR...
		// to insert switch of phenotype from epythelial to mesenchymal according to the value of node EMT

/*
	if ( pCell->phenotype.intracellular->has_variable( "EMT" )){
		if (pCell->phenotype.intracellular->get_boolean_variable_value( "EMT" ))
			custom_detach_cells(pCell);
		else
			custom_cell_attach(pCell);
		} */

	return;
}

/* Calculate repulsion/adhesion between agent and ecm according to its local density */
void add_ecm_interaction(Cell* pC, int index_ecm, int index_voxel )
{
	// Check if there is ECM material in given voxel
	//double dens2 = get_microenvironment()->density_vector(index_voxel)[index_ecm];
	double dens = pC->get_microenvironment()->nearest_density_vector(index_voxel)[index_ecm];
	double ecmrad = sqrt(3.0) * pC->get_microenvironment()->mesh.dx * 0.5;
	// if voxel is "full", density is 1
	dens = std::min( dens, 1.0 ); 
	if ( dens > EPSILON )
	{
		// Distance between agent center and ECM voxel center
		pC->displacement = pC->position - pC->get_container()->underlying_mesh.voxels[index_voxel].center;
		double distance = norm(pC->displacement);
		// Make sure that the distance is not zero
		distance = std::max(distance, EPSILON);
		
		double dd = pC->phenotype.geometry.radius + ecmrad;  
		double dnuc = pC->phenotype.geometry.nuclear_radius + ecmrad;  

		double tmp_r = 0;
		// Cell overlap with ECM node, add a repulsion term
		if ( distance < dd )
		{
			// repulsion stronger if nucleii overlap, see Macklin et al. 2012, 2.3.1
			if ( distance < dnuc )
			{
				double M = 1.0;
				double c = 1.0 - dnuc/dd;
				c *= c;
				c -= M;
				tmp_r = c*distance/dnuc + M;
				pC->custom_data["nucleus_deform"] += (dnuc-distance);
			}
			else
			{
				tmp_r = ( 1 - distance / dd );
				tmp_r *= tmp_r;
			}
			tmp_r *= dens * PhysiCell::parameters.doubles("cell_ecm_repulsion");
		}

		// Cell adherence to ECM through integrins
		double max_interactive_distance = (PhysiCell::parameters.doubles("max_interaction_factor")*pC->phenotype.geometry.radius) + ecmrad;
		if ( distance < max_interactive_distance ) 
		{	
			double temp_a = 1 - distance/max_interactive_distance; 
			temp_a *= temp_a; 
			/* \todo change dens with a maximal density ratio ? */

			pC->custom_data["ecm_contact"] += dens * (max_interactive_distance-distance);
			// temp_a *= dens * ( static_cast<Cell*>(this) )->integrinStrength();

			double temp_integrins = get_integrin_strength( pC->custom_data["pintegrin"] );

			temp_a *= dens * temp_integrins;
			
			tmp_r -= temp_a;
		}
		
		/////////////////////////////////////////////////////////////////
		if(tmp_r==0)
			return;
		tmp_r/=distance;

		pC->velocity += tmp_r * pC->displacement;
	}

}

/** \brief (De)-Activate ECM degradation by the cell */
void set_mmp(Cell* pCell, int activate )
{ 
	int ecm_index = pCell->get_microenvironment()->find_density_index("ecm");

	if (activate){
		pCell->phenotype.secretion.uptake_rates[ecm_index] = 0.05 + (PhysiCell::parameters.doubles("ecm_degradation") * pCell->custom_data["pintegrin"]); // this is  alittle bias so that when mmps is activated the degradation is not totally zero !
	}
		
	else
		pCell->phenotype.secretion.uptake_rates[ecm_index] = 0;
	
}

void set_oxygen_motility(Cell* pC, bool active)
{	
	//phenotype.motility.is_motile = active;
		
	if (active){
		pC->phenotype.motility.chemotaxis_index = pC->get_microenvironment()->find_density_index( "oxygen");
		// bias direction is gradient for the indicated substrate 
		pC->phenotype.motility.migration_bias_direction = pC->nearest_gradient(pC->phenotype.motility.chemotaxis_index);
		pC->phenotype.motility.migration_bias = PhysiCell::parameters.doubles("migration_bias");
		pC->phenotype.motility.chemotaxis_direction = 1.0;
		pC->phenotype.motility.migration_speed = PhysiCell::parameters.doubles("migration_speed");
		pC->phenotype.motility.persistence_time = PhysiCell::parameters.doubles("persistence");
		// move up or down gradient based on this direction 
		pC->phenotype.motility.migration_bias_direction *= pC->phenotype.motility.chemotaxis_direction * get_motility_amplitude(pC->custom_data["pmotility"]); 
	}
	else{
		//restore to default
		pC->phenotype.motility.chemotaxis_index = cell_defaults.phenotype.motility.chemotaxis_index; 
		pC->phenotype.motility.migration_bias_direction = cell_defaults.phenotype.motility.migration_bias_direction;
		pC->phenotype.motility.migration_bias = cell_defaults.phenotype.motility.migration_bias;
		pC->phenotype.motility.chemotaxis_direction = cell_defaults.phenotype.motility.chemotaxis_direction;
		pC->phenotype.motility.migration_speed = cell_defaults.phenotype.motility.migration_speed;
		pC->phenotype.motility.persistence_time = cell_defaults.phenotype.motility.persistence_time;
	 	}
};

void custom_detach_cells(Cell* pCell){

	double junction = pCell->custom_data["padhesion"];

	for (auto cell_attached : pCell->state.attached_cells){

		double otherJunction = cell_attached->custom_data["padhesion"];

		if (junction * otherJunction < PhysiCell::parameters.doubles("cell_junctions_detach_threshold"))
			detach_cells( cell_attached , pCell );

	};

	if (junction == 0.0)
		pCell->remove_all_attached_cells();

}

void custom_cell_attach(Cell* pCell){

	// check the neighbors of current cells, if there is a high level of integrins, it triggers the attachment

	std::vector<Cell*> neigh = pCell->nearby_interacting_cells(); 

	double junction = pCell->custom_data["padhesion"];

	if (neigh.size() == 0)
		return ;

	for (int i = 0; i != neigh.size(); i++){

		double otherJunction = neigh[i]->custom_data["padhesion"];

		//should I also check the distance between the cells?

		if (junction * otherJunction > PhysiCell::parameters.doubles("cell_junctions_attach_threshold"))
			attach_cells( neigh[i] , pCell ); 

	}

}


void custom_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt)
{
	if( pCell->functions.add_cell_basement_membrane_interactions )
	{
		pCell->functions.add_cell_basement_membrane_interactions(pCell, phenotype,dt);
	}
	
	pCell->state.simple_pressure = 0.0; 
	pCell->state.neighbors.clear(); // new 1.8.0
	
	//First check the neighbors in my current voxel
	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for(neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{
		pCell->add_potentials(*neighbor);
	}
	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();

	for( neighbor_voxel_index = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
		neighbor_voxel_index != neighbor_voxel_index_end; 
		++neighbor_voxel_index )
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
			continue;
		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();
		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
		{
			pCell->add_potentials(*neighbor);
		}
	}

	pCell->custom_data["ecm_contact"] = 0.0;
	pCell->custom_data["nucleus_deform"] = 0.0;

	int ecm_index = BioFVM::microenvironment.find_density_index("ecm");
	if ( ecm_index >= 0 ){
		add_ecm_interaction( pCell, ecm_index, pCell->get_current_mechanics_voxel_index() );
		//add_TGFbeta_interaction(pCell, pCell->get_current_mechanics_voxel_index());
	}

	pCell->update_motility_vector(dt); 
	pCell->velocity += phenotype.motility.motility_vector; 
	
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

void SVG_plot_ecm( std::string filename , Microenvironment& M, double z_slice , double time, std::vector<std::string> (*cell_coloring_function)(Cell*), std::string sub )
{


	double X_lower = M.mesh.bounding_box[0];
	double X_upper = M.mesh.bounding_box[3];
 
	double Y_lower = M.mesh.bounding_box[1]; 
	double Y_upper = M.mesh.bounding_box[4]; 

	double plot_width = X_upper - X_lower; 
	double plot_height = Y_upper - Y_lower; 

	double font_size = 0.025 * plot_height; // PhysiCell_SVG_options.font_size; 
	double top_margin = font_size*(.2+1+.2+.9+.5 ); 

	// open the file, write a basic "header"
	std::ofstream os( filename , std::ios::out );
	if( os.fail() )
	{ 
		std::cout << std::endl << "Error: Failed to open " << filename << " for SVG writing." << std::endl << std::endl; 

		std::cout << std::endl << "Error: We're not writing data like we expect. " << std::endl
		<< "Check to make sure your save directory exists. " << std::endl << std::endl
		<< "I'm going to exit with a crash code of -1 now until " << std::endl 
		<< "you fix your directory. Sorry!" << std::endl << std::endl; 
		exit(-1); 
	} 
	
	Write_SVG_start( os, plot_width , plot_height + top_margin );

	// draw the background 
	Write_SVG_rect( os , 0 , 0 , plot_width, plot_height + top_margin , 0.002 * plot_height , "white", "white" );

	// write the simulation time to the top of the plot
 
	char* szString; 
	szString = new char [1024]; 
 
	int total_cell_count = all_cells->size(); 
 
	double temp_time = time; 

	std::string time_label = formatted_minutes_to_DDHHMM( temp_time ); 
 
	sprintf( szString , "Current time: %s, z = %3.2f %s", time_label.c_str(), 
		z_slice , PhysiCell_SVG_options.simulation_space_units.c_str() ); 
	Write_SVG_text( os, szString, font_size*0.5,  font_size*(.2+1), 
		font_size, PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );
	sprintf( szString , "%u agents" , total_cell_count ); 
	Write_SVG_text( os, szString, font_size*0.5,  font_size*(.2+1+.2+.9), 
		0.95*font_size, PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );
	
	delete [] szString; 


	// add an outer "g" for coordinate transforms 
	
	os << " <g id=\"tissue\" " << std::endl 
	   << "    transform=\"translate(0," << plot_height+top_margin << ") scale(1,-1)\">" << std::endl; 
	   
	// prepare to do mesh-based plot (later)
	
	double dx_stroma = M.mesh.dx; 
	double dy_stroma = M.mesh.dy; 
	double dz_stroma = M.mesh.dz;
	
	os << "  <g id=\"ECM\">" << std::endl; 
  
	int ratio = 1; 
	double voxel_size = dx_stroma / (double) ratio ; 
  
	double half_voxel_size = voxel_size / 2.0; 
	double normalizer = 78.539816339744831 / (voxel_size*voxel_size*voxel_size); 
	
	// color the dark background 

	for (int n = 0; n < M.number_of_voxels(); n++)
	{
		auto current_voxel = M.voxels(n);
		int z_center = current_voxel.center[2];
		double z_displ = z_center -  dz_stroma/2;
		
		double z_compare = z_displ;

		if (default_microenvironment_options.simulate_2D == true){
		z_compare = z_center;
		};

		if (z_slice == z_compare){
			int x_center = current_voxel.center[0];
			int y_center = current_voxel.center[1];
			
			double x_displ = x_center -  dx_stroma/2;
			double y_displ = (y_center - dy_stroma) +  dy_stroma/2;
			//std::cout <<  x_displ - X_lower << "  __  " << y_displ - Y_lower << "\n" ;
			int sub_index = M.find_density_index(sub);

			double concentration = M.density_vector(n)[sub_index];

			std::vector< std::string > output( 4 , "black" );
			int color = (int) round( concentration * 255.0 );
			if(color > 255){
				color = 255;
			}
			char szTempString [128];

			int green = 255 - color;
			int blue = 255 - color;
			int red = 255 - color;


			sprintf( szTempString , "rgb(%u, %u, %u)", red, green, blue );
			output[0].assign( szTempString );

			Write_SVG_rect( os , x_displ - X_lower , y_displ - Y_lower, dx_stroma, dy_stroma , 0 , "none", output[0] );
		}

	}

	// Write_SVG_rect( os , 0 , 0 , plot_width, plot_height , 0 , "none", "black" );

 
 // color in the background ECM
/* 
 if( ECM.TellRows() > 0 )
 {
  // find the k corresponding to z_slice
  
  
  
  Vector position; 
  *position(2) = z_slice; 
  

  // 25*pi* 5 microns^2 * length (in source) / voxelsize^3
  
  for( int j=0; j < ratio*ECM.TellCols() ; j++ )
  {
   // *position(1) = *Y_environment(j); 
   *position(1) = *Y_environment(0) - dy_stroma/2.0 + j*voxel_size + half_voxel_size; 
   
   for( int i=0; i < ratio*ECM.TellRows() ; i++ )
   {
    // *position(0) = *X_environment(i); 
    *position(0) = *X_environment(0) - dx_stroma/2.0 + i*voxel_size + half_voxel_size; 
	
    double E = evaluate_Matrix3( ECM, X_environment , Y_environment, Z_environment , position );	
	double BV = normalizer * evaluate_Matrix3( OxygenSourceHD, X_environment , Y_environment, Z_environment , position );
	if( isnan( BV ) )
	{ BV = 0.0; }

	vector<string> Colors;
	Colors = hematoxylin_and_eosin_stroma_coloring( E , BV );
	Write_SVG_rect( os , *position(0)-half_voxel_size-X_lower , *position(1)-half_voxel_size+top_margin-Y_lower, 
	voxel_size , voxel_size , 1 , Colors[0], Colors[0] );
   
   }
  }
 
 }
*/
	os << "  </g>" << std::endl; 

 
	// plot intersecting cells 
	os << "  <g id=\"cells\">" << std::endl; 


	double lowest_pressure = 1.0;
	double max_pressure = 1.0;
	for( int i=0 ; i < total_cell_count ; i++ )
	{
		Cell* pC = (*all_cells)[i]; // global_cell_list[i]; 

		double pressure = pC->custom_data["padhesion"];

		if (pressure < lowest_pressure)
			lowest_pressure = pressure;
		if (pressure > max_pressure)
			max_pressure = pressure;

		static std::vector<std::string> Colors; 
		if( fabs( (pC->position)[2] - z_slice ) < pC->phenotype.geometry.radius )
		{
			double r = pC->phenotype.geometry.radius ; 
			double rn = pC->phenotype.geometry.nuclear_radius ; 
			double z = fabs( (pC->position)[2] - z_slice) ; 
   
			Colors = cell_coloring_function( pC ); 

			os << "   <g id=\"cell" << pC->ID << "\">" << std::endl; 
  
			// figure out how much of the cell intersects with z = 0 
   
			double plot_radius = sqrt( r*r - z*z ); 

			Write_SVG_circle( os, (pC->position)[0]-X_lower, (pC->position)[1]-Y_lower, 
				plot_radius , 0.5, Colors[1], Colors[0] ); 

			// plot the nucleus if it, too intersects z = 0;
			if( fabs(z) < rn && PhysiCell_SVG_options.plot_nuclei == true )
			{   
				plot_radius = sqrt( rn*rn - z*z ); 
			 	Write_SVG_circle( os, (pC->position)[0]-X_lower, (pC->position)[1]-Y_lower, 
					plot_radius, 0.5, Colors[3],Colors[2]); 
			}					  
			os << "   </g>" << std::endl;
		}
	}

	std::cout << max_pressure << "  --  " << lowest_pressure;

	os << "  </g>" << std::endl; 
	
	// end of the <g ID="tissue">
	os << " </g>" << std::endl; 
 
	// draw a scale bar
 
	double bar_margin = 0.025 * plot_height; 
	double bar_height = 0.01 * plot_height; 
	double bar_width = PhysiCell_SVG_options.length_bar; 
	double bar_stroke_width = 0.001 * plot_height; 
	
	std::string bar_units = PhysiCell_SVG_options.simulation_space_units; 
	// convert from micron to mm
	double temp = bar_width;  

	if( temp > 999 && std::strstr( bar_units.c_str() , PhysiCell_SVG_options.mu.c_str() )   )
	{
		temp /= 1000;
		bar_units = "mm";
	}
	// convert from mm to cm 
	if( temp > 9 && std::strcmp( bar_units.c_str() , "mm" ) == 0 )
	{
		temp /= 10; 
		bar_units = "cm";
	}
	
	szString = new char [1024];
	sprintf( szString , "%u %s" , (int) round( temp ) , bar_units.c_str() );
 
	Write_SVG_rect( os , plot_width - bar_margin - bar_width  , plot_height + top_margin - bar_margin - bar_height , 
		bar_width , bar_height , 0.002 * plot_height , "rgb(255,255,255)", "rgb(0,0,0)" );
	Write_SVG_text( os, szString , plot_width - bar_margin - bar_width + 0.25*font_size , 
		plot_height + top_margin - bar_margin - bar_height - 0.25*font_size , 
		font_size , PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() ); 
	
	delete [] szString; 

	// plot runtime 
	szString = new char [1024]; 
	RUNTIME_TOC(); 
	std::string formatted_stopwatch_value = format_stopwatch_value( runtime_stopwatch_value() );
	Write_SVG_text( os, formatted_stopwatch_value.c_str() , bar_margin , top_margin + plot_height - bar_margin , 0.75 * font_size , 
		PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );
	delete [] szString; 

	// draw a box around the plot window
	Write_SVG_rect( os , 0 , top_margin, plot_width, plot_height , 0.002 * plot_height , "rgb(0,0,0)", "none" );
	
	// close the svg tag, close the file
	Write_SVG_end( os ); 
	os.close();
 
	return; 
}
