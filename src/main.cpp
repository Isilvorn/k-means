#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <list>
#include <cmath>
#include <random>

#include "../include/Menu.h"

using namespace std;

default_random_engine              gen;
uniform_real_distribution<double>  distr(0.0,1.0);

#define EPSILON 0.0001

/*
** The Point struct is a container for the x-y pairs that make up the data set.
*/
struct Point {
public:
	Point(void)               { clear();                                } // default constructor
	Point(double a, double b) { data[0] = a;    data[1] = b;    gp = 0; } // alternate constructor #1
	Point(double *c)          { data[0] = c[0]; data[1] = c[1]; gp = 0; } // alternate constructor #2
	~Point(void)              { /* nothing needed */            } // destructor

	void   rand(double=0.0, double=1.0, double=0.0, double=1.0);    // randomize a point
	void   clear(void)        { data[0] = 0.0; data[1] = 0.0; gp = 0.0; }
	double group(void) const  { return gp; }
	void   group(int g)       { gp = g;    }

    Point& operator+=(const Point &rhs) { data[0] += rhs.data[0]; data[1] += rhs.data[1]; }
    Point& operator/=(double d)         { data[0] /= d; data[1] /= d; }
    Point& operator=(const Point &rhs)  { data[0] = rhs.data[0]; data[1] = rhs.data[1]; gp = rhs.gp;     }
    bool   operator==(const Point &rhs) { return ((data[0] == rhs.data[0]) && (data[1] == rhs.data[1])); }
  	double operator[](int n) const      { return (((n>=0) && (n<2))?data[n]:0);                          }

  	friend ostream& operator<<(ostream&,const Point&); // outputs all elements to a stream
  	friend istream& operator>>(istream&, Point&);      // inputs 2 elements from a stream

private:
	double data[2];
	int    gp;
};
// randomize a point
void Point::rand(double xmin, double xmax, double ymin, double ymax) {
	data[0] = xmin + distr(gen) * (xmax-xmin);
	data[1] = ymin + distr(gen) * (ymax-ymin);
}
// input and output points to a stream
ostream& operator<<(ostream& os, const Point& p) 
{ os << p.group() << ": (" << p.data[0] << ", " << p.data[1] << ")"; return os; }
istream& operator>>(istream& is, Point& p) 
{ is >> p.data[0] >> p.data[1]; return is;                       }
// calculate the distance between two points
double dist(Point &p1, Point &p2) { return sqrt((p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1])); }
// compare points using either x or y elements (for use in sorting)
bool compx(const Point &first, const Point &second) { return (first[0] < second[0]); }
bool compy(const Point &first, const Point &second) { return (first[1] < second[1]); }
// get the maximum and minimum x and y values from a list of points
double maxX(list<Point> &plist) { plist.sort(compx); return (plist.back()[0]);  }
double minX(list<Point> &plist) { plist.sort(compx); return (plist.front()[0]); }
double maxY(list<Point> &plist) { plist.sort(compy); return (plist.back()[1]);  }
double minY(list<Point> &plist) { plist.sort(compy); return (plist.front()[1]); }

/*
** The ReadInput() function reads in an input file containing simple x-y pairs.
*/
bool ReadInput(string strFile, list<Point> &plist, bool verbose) {
	ifstream infile;
	Point    pin;

	infile.open(strFile);
	if (infile.is_open()) {
		if (verbose) cout << "Reading file (" << strFile << ")...";
		while (!infile.eof()) {
			infile >> pin;
			plist.push_back(pin);
		}
		if (verbose) cout << "Done." << endl;
		return true;
	}
	else {
		cerr << "Error opening file (" << strFile << ")." << endl;
		return false;
	}
}

/*
** The outdata() function formats the data in a way that can be read
** in by the custom html file.
*/
void outdata(list<Point> &plist, double xmin, double xmax, double ymin, double ymax, int K, 
		double SSE, string fname) {
	ofstream              outfile;   // output file
	list<Point>::iterator it;        // iterator
	Point                 p;

	outfile.open(fname);
	if (outfile.is_open()) {
		// writing the main data file
		it = plist.begin();
		int i=0;
		outfile << "function setdata() {" << endl
			<< "var inputvar =" << endl
			<< "[ " << endl
			<< fixed << setprecision(3);
		while (it != plist.end()) { 
			p = *it;
			outfile << "  [ " 
				<< setw(5) << p[0]      << ", "
				<< setw(5) << p[1]      << ", " 
				<< setw(5) << p.group() << " ]";
			it++;
			if (it == plist.end()) outfile << " ];" << endl;
			else                   outfile << ", "  << endl;
		} // end while (it)
		outfile << "return inputvar;" << endl;
		outfile << "}" << endl << endl;

		// writing formatting parameters
		outfile << fixed << setprecision(3)
			<< "function parms() {" << endl
			<< "var inputvar =" << endl
			<< "[ " << xmin << ", " << xmax << ", " << ymin << ", " << ymax << ", " 
			<< K << ", " << SSE << " ];" << endl
			<< "return inputvar;" << endl
			<< "}" << endl;

		outfile.close();
  	} // end if (outfile)
  	else cerr << "Error opening file for output." << endl;
} // end outdata()

/*
** The SSEcalc() function calculates the sum of squared errors using the calculated centroids.
*/
double SSEcalc(list<Point> &plist, list<Point> &clist) {
	list<Point>::iterator it, it2;
	Point p1, p2;
	double d;
	double SSE = 0.0;

	it2 = clist.begin();
	while (it2 != clist.end()) {
		p2 = *it2;
		it = plist.begin();
		while (it != plist.end()) {
			p1 = *it;
			if (p1.group() == p2.group()) {
				d = dist(p1,p2);
				SSE += d*d;
			}
			it++;
		} // end while (it)
		it2++;
	} // end while (it2)
	return SSE;
} // end SSEcalc()

int main(void) {
	Menu   mainmenu;
	int    idx = 0;
	int    K = 3;
	int    iter = 10;
	int    i, c, g;
	bool   flag;
	double xmax, xmin, ymax, ymin, d, dmin;
	string inputFile = "data/A.txt";
	list<Point> pointlist;
	list<Point> centroids;
	list<Point>::iterator it, it2;
	Point p1,p2;

	mainmenu.load("mainmenu.txt");

	do {
		switch(idx) {
			case 1:   // Read the specified input file
				ReadInput(inputFile, pointlist, true);
				break;
			case 2:   // Select the input file called A.txt
				inputFile = "data/A.txt";
				break;
			case 3:   // Select the input file called B.txt
				inputFile = "data/B.txt";
				break;
			case 4:   // Enter an arbritrary input file
				cout << "Enter new file (with path): ";
				cin >> inputFile;
				break;
			case 5:   // Select a new number of clusters to iterate on
				cout << "Enter in a new K (integer): ";
				cin  >> K;
				break;
			case 6:   // Select a new number of maximum iterations
				cout << "Enter in a new max iterations (integer): ";
				cin  >> iter;
				break;
			case 7:   // Iterate to a solution using k-means
				if (pointlist.size() == 0) flag = ReadInput(inputFile, pointlist, false); else flag = true;
				if (flag) {
					// generating a random set of centroids
					centroids.erase(centroids.begin(),centroids.end());
					xmax = maxX(pointlist);
					xmin = minX(pointlist);
					ymax = maxY(pointlist);
					ymin = minY(pointlist);
					for (int i=0; i < K; i++) {
						p1.rand(xmin,xmax,ymin,ymax);
						it = pointlist.begin();
						dmin = sqrt((xmax-xmin)*(xmax-xmin) + (ymax-ymin)*(ymax-ymin));
						// picking the closest point in the data set to the random point
						// to help avoid empty clusters
						while (it != pointlist.end()) {
							d = dist(p1,(*it));
							if (d < dmin) { dmin = d; p2 = *it; p2.group(i+1); }
							it++;
						}
						centroids.push_back(p2);
					} // end for (i)

					// iterate over max iterations to a solution
					for (i=0; i < iter; i++) {
						// updating data point groups according to which centroid is closest
						// to the point
						it = pointlist.begin();
						while (it != pointlist.end()) {
							p1 = *it;
							it2 = centroids.begin();
							dmin = sqrt((xmax-xmin)*(xmax-xmin) + (ymax-ymin)*(ymax-ymin));
							while (it2 != centroids.end()) {
								p2 = *it2;
								d = dist(p1,p2);
								if (d < dmin) { dmin = d; g = p2.group(); }
								it2++;
							} // end while (it2)
							(*it).group(g);
							it++;
						} // end while (it)

						// updating centroids according to the center of gravity of the groups
						flag = true;
						it2  = centroids.begin();
						while (it2 != centroids.end()) {
							p2 = *it2;
							p1.clear();
							p1.group(p2.group());
							c = 0;
							it = pointlist.begin();
							while (it != pointlist.end()) {
								if ((*it).group() == p2.group()) {
									c++;
									p1 += *it;
								} // end if (it)
								it++;
							} // end while (it)
							p1 /= c;
							if (dist((*it2),p1) > EPSILON) { (*it2) = p1; flag=false; }
							it2++;
						} // end while (it2)

						// break out of the loop if all of the centroids are the same from the
						// last iteration
						if (flag) break; 
					} // end for (i) 
					cout << "Converged after " << i << " iterations." << endl;
				} // end if (flag)
				break;
			case 8:   // Assemble agglomerative heirarchichal clusters
				break;
			case 9:   // Display the data points by sending them to stdout
				cout << "Data Points" << endl;
				cout << "===========" << endl;
				it = pointlist.begin();
				while (it != pointlist.end()) {
					cout << *it << endl;
					it++;
				}
				break;
			case 10:   // Display the centroids by sending them to stdout
				cout << "Centroids" << endl;
				cout << "=========" << endl;
				it = centroids.begin();
				while (it != centroids.end()) {
					cout << *it << endl;
					it++;
				}
				cout << endl << "SSE = " << SSEcalc(pointlist, centroids) << endl << endl;				
				break;
			case 11:  // write data to file
				if (pointlist.size() > 0) {
					xmax = maxX(pointlist);
					xmin = minX(pointlist);
					ymax = maxY(pointlist);
					ymin = minY(pointlist);
					outdata(pointlist, xmin, xmax, ymin, ymax, K, SSEcalc(pointlist, centroids), "setdata.js");
				} // end if (pointlist.size())
				else cerr << "Nothing to write." << endl;
				break;
		}

		mainmenu.draw(0,40);
		mainmenu.addenda("K = ",K,2,false);
		mainmenu.addenda("Max iterations = ",iter,5,false);
		mainmenu.addenda("Current File: " + inputFile,true);

	 } while ((idx=mainmenu.prompt()) != 0); {

	}


}