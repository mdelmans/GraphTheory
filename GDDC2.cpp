/*
GDDC2 - Graphlet Degree Distribution Comparison 2

A tool for finding GDD simmilarity between two groups of cells in a segmentation.

Running:

./GDDC2 segmentation.xml

Author:
M.Delmans

*/

#include <CGAL/Linear_cell_complex.h>
#include <CMPlant/2dCellAttribute.hpp>

typedef CGAL::Linear_cell_complex_traits<2, CGAL::Exact_predicates_inexact_constructions_kernel> Traits;
typedef CGAL::Linear_cell_complex<2, 2, Traits, PlantItem> PCC;

#include <armadillo>
#include <CMPlant/IO.hpp>
#include <CMPlant/Network.hpp> 
#include <CMPlant/2DDraw.hpp>

#include <iterator>
#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/SoQtRenderArea.h>
#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>


std::pair <double, double> gddSimilarity(arma::mat gx1, arma::mat gx2){
	
	int gxMax = int (std::max( gx1.max(), gx2.max() )); 

	arma::umat dgj1 = arma::hist(gx1, arma::linspace<arma::vec>(0, gxMax, gxMax+1));
	arma::umat dgj2 = arma::hist(gx2, arma::linspace<arma::vec>(0, gxMax, gxMax+1));

	std::cout << gx1.n_rows << ":" << gx1.n_cols << "\n";
	std::cout << gx2.n_rows << ":" << gx2.n_cols << "\n"; 

	std::cout << dgj1.n_rows << ":" << dgj1.n_cols << "\n"; 
	std::cout << dgj2.n_rows << ":" << dgj2.n_cols << "\n"; 


	for (int col = 0; col < dgj1.n_cols; col ++){
		for (int row = 1; row < dgj1.n_cols; row ++){
			dgj1(row, col) /= row;
		}
	}

	for (int col = 0; col < dgj2.n_cols; col ++){
		for (int row = 1; row < dgj2.n_cols; row ++){
			dgj2(row, col) /= row;
		}
	}

	dgj1(0, arma::span::all) = arma::zeros<arma::umat>(1, dgj1.n_cols);
	dgj2(0, arma::span::all) = arma::zeros<arma::umat>(1, dgj2.n_cols);

	arma::urowvec tgj1 = arma::sum(dgj1);
	arma::urowvec tgj2 = arma::sum(dgj2);

	std::cout << tgj1.n_elem << " " << tgj2.n_elem << "\n";

	for (int col = 0; col < dgj1.n_cols; col ++){
		
		if (tgj1(col) != 0){
			dgj1.col(col) /= tgj1(col);
		}
		else{
			dgj1.col(col) = arma::zeros<arma::uvec>( dgj1.n_rows );
		}

		if (tgj2(col) != 0){
			dgj2.col(col) /= tgj2(col);
		}
		else{
			dgj2.col(col) = arma::zeros<arma::uvec>( dgj2.n_rows );
		}
	}

	arma::vec D = arma::zeros<arma::vec>(dgj1.n_cols);

	for(int i = 0; i < dgj1.n_cols; i++){
		D(i) += sqrt( arma::sum( arma::pow(dgj1.col(i) - dgj2.col(i), 2) ) ) / 1.4142;
	}

	arma::vec A = arma::ones<arma::vec>(dgj1.n_cols) - D;

	double mu = arma::mean(A);
	double geoMu = pow(arma::prod(A), 1.0/73);

	return {mu, geoMu};
}	


int main(int argc, char **argv){

	PCC pcc;
	CMPlant::importFromSegmentation(pcc, argv[1], 1);
	std::cout << "Imported " << argv[1] << "\n";
	runNCount(pcc);
	std::cout << "Counted graphlet signatures" << "\n";

	int nPoints(0);
	int nVars(73);

	arma::mat gx, gx1, gx2;
	std::unordered_map<int, std::vector<double>> gxMap;
	int nCells = 0;

	forEachCell(pcc, [&](PCC& pcc, CellIterator it){
		std::vector<double> dDist( pcc.info<2>(it).graphletDist.begin(), pcc.info<2>(it).graphletDist.end() );

		gxMap[pcc.info<2>(it).id] = dDist;
		nCells++;
		if (pcc.info<2>(it).group == 1){
			gx1.insert_rows(gx1.n_rows, arma::rowvec(dDist));
		}
		if (pcc.info<2>(it).group == 2){
			gx2.insert_rows(gx2.n_rows, arma::rowvec(dDist));
		}
	});

	std::pair<double,double> similarity = gddSimilarity(gx1, gx2);

	std::cout << "Aarithm=" << similarity.first << "\n";
	std::cout << "Ageo=" << similarity.second << "\n\n";


	arma::rowvec weight({1-log(1.0)/log(73), 1-log(2.0)/log(73), 1-log(2.0)/log(73), 1-log(2.0)/log(73), 1-log(3.0)/log(73), 1-log(4.0)/log(73), 1-log(3.0)/log(73), 1-log(3.0)/log(73), 1-log(4.0)/log(73), 1-log(3.0)/log(73), 1-log(4.0)/log(73), 1-log(4.0)/log(73), 1-log(4.0)/log(73), 1-log(4.0)/log(73), 1-log(3.0)/log(73), 1-log(4.0)/log(73), 1-log(6.0)/log(73), 1-log(5.0)/log(73), 1-log(4.0)/log(73), 1-log(5.0)/log(73), 1-log(6.0)/log(73), 1-log(6.0)/log(73), 1-log(4.0)/log(73), 1-log(4.0)/log(73), 1-log(4.0)/log(73), 1-log(5.0)/log(73), 1-log(7.0)/log(73), 1-log(4.0)/log(73), 1-log(6.0)/log(73), 1-log(6.0)/log(73), 1-log(7.0)/log(73), 1-log(4.0)/log(73), 1-log(6.0)/log(73), 1-log(6.0)/log(73), 1-log(6.0)/log(73), 1-log(5.0)/log(73), 1-log(6.0)/log(73), 1-log(7.0)/log(73), 1-log(7.0)/log(73), 1-log(5.0)/log(73), 1-log(7.0)/log(73), 1-log(6.0)/log(73), 1-log(7.0)/log(73), 1-log(6.0)/log(73), 1-log(5.0)/log(73), 1-log(5.0)/log(73), 1-log(6.0)/log(73), 1-log(8.0)/log(73), 1-log(7.0)/log(73), 1-log(6.0)/log(73), 1-log(6.0)/log(73), 1-log(8.0)/log(73), 1-log(6.0)/log(73), 1-log(9.0)/log(73), 1-log(5.0)/log(73), 1-log(6.0)/log(73), 1-log(4.0)/log(73), 1-log(6.0)/log(73), 1-log(6.0)/log(73), 1-log(7.0)/log(73), 1-log(8.0)/log(73), 1-log(6.0)/log(73), 1-log(6.0)/log(73), 1-log(8.0)/log(73), 1-log(7.0)/log(73), 1-log(6.0)/log(73), 1-log(7.0)/log(73), 1-log(7.0)/log(73), 1-log(8.0)/log(73), 1-log(5.0)/log(73), 1-log(6.0)/log(73), 1-log(6.0)/log(73), 1-log(4.0) / log(73)});

	for(int i = 0; i<nCells; i++){
		gx.insert_rows(gx.n_rows, arma::rowvec(gxMap[i]) % weight);
	}

	QWidget *window = SoQt::init(argv[0]);
  	if (window == NULL) exit(1);

  	arma::mat score;
  	arma::mat coeff;
  	arma::vec latent;

	princomp(coeff, score, latent, gx);

	std::cout << "Variation coverage:\n";
	latent = latent / arma::sum(latent);
	std::cout << arma::sum( latent.subvec(0,2) ) << "\n\n";

	arma::colvec r = score.col(0);
	arma::colvec g = score.col(1);
	arma::colvec b = score.col(2);

	r = (r - min(r) ) / (max(r) - min(r));
	g = (g - min(g) ) / (max(g) - min(g));
	b = (b - min(b) ) / (max(b) - min(b));

	int cid;

	forEachCell(pcc, [&](PCC& pcc, CellIterator it){
		cid = pcc.info<2>(it).id;
		pcc.info<2>(it).color = SbColor( r(cid), g(cid), b(cid) );
		copyAttributeToCell(pcc, it);
	});


	SoSeparator *root = new SoSeparator;
  	root->ref();

  	SoPerspectiveCamera *camera = new SoPerspectiveCamera;
  	root->addChild(camera);

  	root->addChild(new SoDirectionalLight);

  	SoSeparator* plantScene = makePlantScene(pcc);
  	updateScene(pcc, plantScene);

  	root->addChild(plantScene);

  	SoQtRenderArea *renderArea = new SoQtRenderArea(window);

  	camera->viewAll(root, renderArea->getViewportRegion());

  	renderArea->setSceneGraph(root);
  	renderArea->setTitle("Engine Spin");
  	renderArea->show();

  	SoQt::show(window);
  	SoQt::mainLoop();

  	delete renderArea;
  	root->unref();

  	return 0;
}
