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
{ os << (int)p.group() << ": (" << p.data[0] << ", " << p.data[1] << ")"; return os; }
istream& operator>>(istream& is, Point& p) 
{ is >> p.data[0] >> p.data[1]; return is; }
// calculate the distance between two points
double dist(Point &p1, Point &p2) { return sqrt((p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1])); }
// compare points using x/y elements or the group (for use in sorting)
bool compx(const Point &first, const Point &second) { return (first[0] < second[0]);           }
bool compy(const Point &first, const Point &second) { return (first[1] < second[1]);           }
bool compg(const Point &first, const Point &second) { return (first.group() < second.group()); }
// get the maximum and minimum x and y values from a list of points
double maxX(list<Point> &plist) { plist.sort(compx); return (plist.back()[0]);  }
double minX(list<Point> &plist) { plist.sort(compx); return (plist.front()[0]); }
double maxY(list<Point> &plist) { plist.sort(compy); return (plist.back()[1]);  }
double minY(list<Point> &plist) { plist.sort(compy); return (plist.front()[1]); }

typedef list<Point> ListOfPoints;

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
		double SSE, int mode, string fname) {
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
			<< K << ", " << SSE << ", " << mode << " ];" << endl
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

/*
** The updateCentroids() function takes the pointlist and updates the centroids according
** to the new center of gravity of the entire group.  It returns true if the update resulted
** in a change.  If there was no change it returns false.
*/
bool updateCentroids(list<Point> &plist, list<Point> &clist) {
	bool                  flag = false;
	int                   c;
	list<Point>::iterator it, it2;
	Point                 p1, p2;

	it2  = clist.begin();
	while (it2 != clist.end()) {
		p2 = *it2;
		p1.clear();
		p1.group(p2.group());
		c = 0;
		it = plist.begin();
		while (it != plist.end()) {
			if ((*it).group() == p2.group()) {
				c++;
				p1 += *it;
			} // end if (it)
			it++;
		} // end while (it)
		p1 /= c;
		if (dist((*it2),p1) > EPSILON) { (*it2) = p1; flag=true; }
		it2++;
	} // end while (it2)

	return flag;
} // end updateCentroids()

/*
** The updatePointGroups() function updates the group of each point according to the
** the nearest centroid.  Returns true if a point was updated.  If there was no change
** it returns false.
*/
bool updatePointGroups(list<Point> &plist, list<Point> &clist, double dmax) {
	double                g, d, dmin;
	bool                  flag = false;
	list<Point>::iterator it, it2;
	Point                 p1, p2;

	it = plist.begin();
	while (it != plist.end()) {
		p1 = *it;
		it2 = clist.begin();
		dmin = dmax; // setting to the maximum possible value to start
		while (it2 != clist.end()) {
			p2 = *it2;
			d = dist(p1,p2);
			if (d < dmin) { dmin = d; g = p2.group(); }
			it2++;
		} // end while (it2)
		if ((*it).group() != g) { (*it).group(g); flag = true; }
		it++;
	} // end while (it)

	return flag;
} // end updatePointGroups()

/*
** The randomizeCentroids() function creates a given number of centroids at random places
** within the bounds of the data set.  Returns the length of the diagonal across the
** breadth of the rectangle bounding all data points as an indication of the maximum
** possible distance from point to point in the dataset.
*/
double randomizeCentroids(list<Point> &plist, list<Point> &clist, int K) {
	Point                 p1, p2;
	double                xmin, xmax, ymin, ymax, d, dmin, dmax;
	list<Point>::iterator it;

	clist.erase(clist.begin(),clist.end()); // removing any residual data

	// finding the bounding rectangle of all points in the list
	xmax = maxX(plist);
	xmin = minX(plist);
	ymax = maxY(plist);
	ymin = minY(plist);
	dmax = sqrt((xmax-xmin)*(xmax-xmin) + (ymax-ymin)*(ymax-ymin));


	for (int i=0; i < K; i++) {
		p1.rand(xmin,xmax,ymin,ymax);
		it = plist.begin();
		dmin = dmax; // setting to the maximum possible value to start
		// picking the closest point in the data set to the random point
		// to help avoid empty clusters
		while (it != plist.end()) {
			d = dist(p1,(*it));
			if (d < dmin) { dmin = d; p2 = *it; p2.group(i+1); }
			it++;
		}
		clist.push_back(p2);
	} // end for (i)

	return dmax;
}

/*
** The initCentroids() function initializes all centroids (groups) as singleton clusters.
*/
double initCentroids(list<Point> &plist, list<Point> &clist, list<ListOfPoints> hlist) {
	int                   i;
	double                xmin, xmax, ymin, ymax, dmax;
	list<Point>          *lp;
	list<Point>::iterator it;
	Point                 p1, p2;

	lp = new ListOfPoints;         // allocating a new list
	i = 0;
	clist.erase(clist.begin(),clist.end());
	it = plist.begin();
	while (it != plist.end()) {
		p1 = *it;                  // selecting the next point in the point list
		p1.group(i+1);             // changing the group to an incremented index
		(*lp).push_front(*it);     // pushing it onto the list
		clist.push_front(p1);      // creating a centroid on that point
		(*it).group(i+1);          // recording the group change in the point list
		it++;
		i++;
	}
	hlist.push_front(*lp);         // pushing the new list onto the hclusters list of lists

	// finding the bounding rectangle of all points in the list
	xmax = maxX(plist);
	xmin = minX(plist);
	ymax = maxY(plist);
	ymin = minY(plist);
	dmax = sqrt((xmax-xmin)*(xmax-xmin) + (ymax-ymin)*(ymax-ymin));

	return dmax;

}

/*
** The findMIN function finds the MIN similarity between the group for the supplied point 
** among the other groups.  Return value is the point in the other groups that is the
** closest.  "p1" is the point in "this" group, that is the closest to the point in another
** group, while "p2" is the point in the "other" group that is the closest to "this" group.
*/
void findMIN(list<Point> &plist, Point &p1, Point &p2, double dmax) {
	double                d, dmin;
	int                   g;
	list<Point>::iterator it, it2;
	Point pa, pb;

	it = plist.begin();
	dmin = dmax;
	g    = p1.group();
	// Searching all pairs between the group identified by the supplied group and all other
	// groups to find the pair with the minimum distance.  The point in the other group
	// that meets that criteria is the return value.
	while (it != plist.end()) {
		it2 = plist.begin();
		if ((*it).group() == g) {
			while (it2 != plist.end()) {
				if ((*it2).group() != g) {
					d = dist((*it),(*it2));
					if (d < dmin) {
						dmin = d;
						p1 = *it;
						p2 = *it2;
					} // end if (d)
				} // end if (it2)
				it2++;
			} // end while (it2)
		} // end if (it)
		it++;
	} // end while (it)

} // end findMIN()

/*
** The findMAX function finds the MAX similarity between the group for the supplied point 
** among the other groups.  Return value is the point in the other groups that is the
** smallest of the furthest distances between points in the groups.  "p1" is the point in "this" 
** group, that is the furthest to the point in the "other" group, while "p2" is the point in the 
** "other" group that is the furthest to the point in "this" group.
*/
void findMAX(list<Point> &plist, list<Point> &clist, Point &p1, Point &p2, double dmax) {
	double                d, dmin, dmaxg;
	int                   g;
	bool                  eg, eng;
	list<Point>::iterator it, it2, it3;
	Point pa, pb;

	dmin = dmax;
	g    = p1.group();
	plist.sort(compg);  // sorting by group to help speed up the iterating process

	// Searching all pairs between the group identified by the supplied group and all other
	// groups to find the pair with the minimum of the maximum-within-group distance.  The point 
	// in the other group that meets that criteria is the return value.
	it = clist.begin();
	while (it != clist.end()) {

		if ((*it).group() != g) {
			eg = false;
			it2 = plist.begin();
			while (it2 != plist.end()) {
				if ((*it2).group() == g) {
					eg = true;
					dmaxg = 0.0;
					eng = false;
					it3 = plist.begin();
					while (it3 != plist.end()) {
						if ((*it3).group() == (*it).group()) {
							eng = true;
							d = dist((*it2),(*it3));
							if (d > dmaxg) {
								dmaxg = d;
								pa = *it2;
								pb = *it3;
							} // end if (d)
						} // end if (it3)
						else if (eng) break;
						it3++;
					} // end while (it3)
				} // end if (it2)
				else if (eg) break;
				it2++;
			} // end while (it2)
			// checking the distance between the furthest points between group "g" and 
			// group (*it).group()
			d = dist(pa, pb);
			if (d < dmin) {
				dmin = d;
				p1 = pa;
				p2 = pb;
			} // end if (d)
		} // end if (it)
		it++;
	} // end while (it)
} // end findMIN()

/*
** The findGAVG function finds the group average similarity between the group for the supplied point 
** among the other groups.  Return value is a point in the other groups that is the
** smallest average distance between points in the groups.  "p1" is the point in "this" 
** group, that is the furthest to the point in the "other" group, while "p2" is the point in the 
** "other" group that is the furthest to the point in "this" group.
*/
void findGAVG(list<Point> &plist, list<Point> &clist, Point &p1, Point &p2, double &dmin, double dmax) {
	double                d, davg;
	int                   g, n;
	bool                  eg, eng;
	list<Point>::iterator it, it2, it3;
	Point pa, pb;

	dmin = dmax;
	g    = p1.group();
	plist.sort(compg);  // sorting by group to help speed up the iterating process

	// Searching all pairs between the group identified by the supplied group and all other
	// groups to find the pair with the minimum average-between-group distance.  The point 
	// in the other group that meets that criteria is the return value.
	it = clist.begin();
	while (it != clist.end()) {

		if ((*it).group() != g) {
			eg   = false;
			davg = 0.0;
			n    = 0;
			it2  = plist.begin();
			while (it2 != plist.end()) {
				if ((*it2).group() == g) {
					eg   = true;
					eng  = false;
					it3  = plist.begin();
					while (it3 != plist.end()) {
						if ((*it3).group() == (*it).group()) {
							eng = true;
							d = dist((*it2),(*it3));
							davg += d;
							n++;
						} // end if (it3)
						else if (eng) break;
						it3++;
					} // end while (it3)
				} // end if (it2)
				else if (eg) break;
				it2++;
			} // end while (it2)
			// checking the distance between the average distance between group "g" and 
			// group (*it).group()
			davg /= n;
			if (davg < dmin) {
				dmin = davg;
				p2.group((*it).group());
			} // end if (d)
		} // end if (it)
		it++;
	} // end while (it)
} // end findGAVG()

/*
** The findCENTR function finds the similarity between the group for the supplied point 
** among the other groups by comparing their centroids.  Return value is the centroid of the other group 
** that has the smallest distance between group centroids.  "p1" is the centroid of "this" 
** group, while "p2" is the centroid of the "other" group.
*/
void findCENTR(list<Point> &plist, list<Point> &clist, Point &p1, Point &p2, double dmax) {
	double                d, dmin;
	int                   g, n;
	bool                  eg, eng;
	list<Point>::iterator it, it2, it3;
	Point pa, pb;

	dmin = dmax;
	g    = p1.group();

	// ensuring all centroids are up-to-date
	updateCentroids(plist, clist);
	// finding the centroid corresponding to the group of the point "p1"
	it = clist.begin();
	while (it != clist.end()) {
		if ((*it).group() == p1.group()) {
			p1 = *it;
			break;
		}
		it++;
	}

	// Searching all centroid pairs between the group identified by "p1" and all other
	// groups to find the pair with the minimum centroid-to-centroid distance.  The point 
	// in the other group that meets that criteria is the return value.
	it = clist.begin();
	while (it != clist.end()) {

		if ((*it).group() != g) {
			d = dist(p1, *it);
			if (d < dmin) {
				dmin = d;
				p2 = *it;
			} // end if (d)
		} // end if (it)
		it++;

	} // end while (it)
} // end findCENTR()

/*
** The mergeGroups() function takes the pointlist and two groups numbers and inputs, and then
** merges the groups together.  The resulting group number is the lowest magnitude group number.
*/
void mergeGroups(list<Point> &plist, list<Point> &clist, int g1, int g2) {
	list<Point>::iterator it;
	int g;

	// "g" is the lowest magnitude group number, and is to be the newly merged group number
	g = (g1 <= g2)?g1:g2;

	// setting all of the points in the numerically larger group to the new group number
	it = plist.begin();
	while (it != plist.end()) {
		if (((*it).group() == g1) || ((*it).group() == g2)) {
			if ((*it).group() != g) (*it).group(g);
		}
		it++;
	}

	// "g" is now the group to remove
	g = (g1 > g2)?g1:g2;

	// removing the old group number from the list of groups
	it = clist.begin();
	while (it != clist.end()) {
		if ((*it).group() == g) {
			it = clist.erase(it);
			break;
		}
		it++;
	}
}

/*
** The minMerge() function finds the two groups with the nearest points to each other and merges
** them.
*/ 
void minMerge(list<Point> &plist, list<Point> &clist, list<ListOfPoints> &hlist, double dmax) {
	double                 d, dmin;
	list<Point>::iterator  it;
	list<Point>           *lp;
	Point                  p1, p2, pa, pb;

	// find out which two groups are closest
	dmin = dmax;
	it = clist.begin();
	while (it != clist.end()) {
		p1 = *it;
		findMIN(plist, p1, p2, dmax);
		d = dist(p1, p2);
		if (d < dmin) {
			dmin = d;
			pa = p1;
			pb = p2;
		} // end if (d)
		it++;
	} // end while (it)

	// merge the two groups
	mergeGroups(plist, clist, pa.group(), pb.group());

	// add the new grouping to the heirarchical cluster list
	lp = new ListOfPoints;     // allocating a new list
	it = plist.begin();
	while (it != plist.end()) {
		(*lp).push_front(*it); // pushing it onto the list
		it++;
	} // end while (it)
	hlist.push_front(*lp);     // pushing the new list onto the hclusters list of lists


} // end minMerge()

/*
** The maxMerge() function finds the two groups with the smallest value of distance between 
** points that are the maximum distance from each other within the group, and then merges
** them.
*/ 
void maxMerge(list<Point> &plist, list<Point> &clist, list<ListOfPoints> &hlist, double dmax) {
	double                 d, dmin;
	list<Point>::iterator  it;
	list<Point>           *lp;
	Point                  p1, p2, pa, pb;

	// find out which two groups are closest
	dmin = dmax;
	it = clist.begin();
	while (it != clist.end()) {
		p1 = *it;
		findMAX(plist, clist, p1, p2, dmax);
		d = dist(p1, p2);
		if (d < dmin) {
			dmin = d;
			pa = p1;
			pb = p2;
		} // end if (d)
		it++;
	} // end while (it)

	// merge the two groups
	mergeGroups(plist, clist, pa.group(), pb.group());

	// add the new grouping to the heirarchical cluster list
	lp = new ListOfPoints;     // allocating a new list
	it = plist.begin();
	while (it != plist.end()) {
		(*lp).push_front(*it); // pushing it onto the list
		it++;
	} // end while (it)
	hlist.push_front(*lp);     // pushing the new list onto the hclusters list of lists

} // end maxMerge()

/*
** The gavgMerge() function finds the two groups with the smallest difference between their
** point-to-point group averages.
*/ 
void gavgMerge(list<Point> &plist, list<Point> &clist, list<ListOfPoints> &hlist, double dmax) {
	double                 d, dmin;
	list<Point>::iterator  it;
	list<Point>           *lp;
	Point                  p1, p2, pa, pb;

	// find out which two groups are closest
	dmin = dmax;
	it = clist.begin();
	while (it != clist.end()) {
		p1 = *it;
		findGAVG(plist, clist, p1, p2, d, dmax);
		if (d < dmin) {
			dmin = d;
			pa = p1;
			pb = p2;
		} // end if (d)
		it++;
	} // end while (it)

	// merge the two groups
	mergeGroups(plist, clist, pa.group(), pb.group());

	// add the new grouping to the heirarchical cluster list
	lp = new ListOfPoints;     // allocating a new list
	it = plist.begin();
	while (it != plist.end()) {
		(*lp).push_front(*it); // pushing it onto the list
		it++;
	} // end while (it)
	hlist.push_front(*lp);     // pushing the new list onto the hclusters list of lists

} // end gavgMerge()

/*
** The centrMerge() function finds the two groups with the smallest difference between their
** centroids.
*/ 
void centrMerge(list<Point> &plist, list<Point> &clist, list<ListOfPoints> &hlist, double dmax) {
	double                 d, dmin;
	list<Point>::iterator  it;
	list<Point>           *lp;
	Point                  p1, p2, pa, pb;

	// find out which two groups are closest
	dmin = dmax;
	it = clist.begin();
	while (it != clist.end()) {
		p1 = *it;
		findCENTR(plist, clist, p1, p2, dmax);
		d = dist(p1, p2);
		if (d < dmin) {
			dmin = d;
			pa = p1;
			pb = p2;
		} // end if (d)
		it++;
	} // end while (it)

	// merge the two groups
	mergeGroups(plist, clist, pa.group(), pb.group());

	// add the new grouping to the heirarchical cluster list
	lp = new ListOfPoints;     // allocating a new list
	it = plist.begin();
	while (it != plist.end()) {
		(*lp).push_front(*it);     // pushing it onto the list
		it++;
	} // end while (it)
	hlist.push_front(*lp);     // pushing the new list onto the hclusters list of lists

} // end gavgMerge()

int main(void) {
	Menu   mainmenu;
	int    idx = 0;
	int    K = 3;
	int    iter = 10;
	int    i, c, g, mode = 0;
	bool   flag;
	double xmax, xmin, ymax, ymin, d, dmin, dmax;
	string inputFile = "data/A.txt";
	ListOfPoints       pointlist;  // points read in from data file
	ListOfPoints       centroids;  // centroids of all clusters (could also be considerd a list of the groups)
	ListOfPoints       *lp;        // temp lists for use in heirarchical clusters
	list<ListOfPoints> hclusters;  // heirarchical cluster lists

	list<Point>::iterator        it, it2;
	list<ListOfPoints>::iterator lit1, lit2;
	Point p1, p2, pa, pb;

	mainmenu.load("mainmenu.txt");

	do {
		switch(idx) {
			case 1:   // Read the specified input file
				pointlist.erase(pointlist.begin(), pointlist.end());
				hclusters.erase(hclusters.begin(), hclusters.end());
				centroids.erase(centroids.begin(), centroids.end());
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

				// checking to see if a file has already been read in, if not it reads in the selected file
				if (pointlist.size() == 0) flag = ReadInput(inputFile, pointlist, false); else flag = true;

				if (flag) {
					// generating a random set of centroids, dmax is the length of the diagonal of the
					// bounding rectangle for all points in the data set (which gives an upper bound
					// of the distance that can be expected between points)
					dmax = randomizeCentroids(pointlist, centroids, K);

					// iterate over max iterations to a solution
					for (i=0; i < iter; i++) {
						// updating data point groups according to which centroid is closest
						// to the point
						flag = updatePointGroups(pointlist, centroids, dmax);

						if (flag) { // only bothers updating centroids if any points changed clusters
							// updating centroids according to the center of gravity of the groups
							flag = updateCentroids(pointlist, centroids);
						} // end if (flag)

						// break out of the loop if all of the centroids are the same as from the
						// last iteration
						if (!flag) break; 
					} // end for (i) 
					cout << "  Converged after " << i << " iterations." << endl;
					cout << setprecision(3) << fixed;
					cout << "  SSE = " << SSEcalc(pointlist, centroids) << endl;
					mode = 0;

				} // end if (flag)
				break;
			case 8:   // Assemble agglomerative heirarchichal clusters using MIN

				// checking to see if a file has already been read in, if not it reads in the selected file
				if (pointlist.size() == 0) flag = ReadInput(inputFile, pointlist, false); else flag = true;

				if (flag) {
					// assembling the initial singleton groups (each point is its own group)
					dmax = initCentroids(pointlist, centroids, hclusters);

					// iterate through until the desired number of clusters is reached (K)
					for (int i=0; i < (pointlist.size() - K); i++)
						minMerge(pointlist, centroids, hclusters, dmax);

					// calculate the centroids of each group so that the SSE calculation is correct
					updateCentroids(pointlist, centroids);

					// sort lists by the group number for aesthetic reasons
					pointlist.sort(compg);
					centroids.sort(compg);

					cout << "  Completed after " << (int)(pointlist.size()-K) << " merge steps." << endl;
					cout << setprecision(3) << fixed;
					cout << "  SSE = " << SSEcalc(pointlist, centroids) << endl;
					mode = 1;

				} // end if (flag)

				break;
			case 9:   // Assemble agglomerative heirarchichal clusters using MAX

				// checking to see if a file has already been read in, if not it reads in the selected file
				if (pointlist.size() == 0) flag = ReadInput(inputFile, pointlist, false); else flag = true;

				if (flag) {
					pointlist.sort(compy);

					// assembling the initial singleton groups (each point is its own group)
					dmax = initCentroids(pointlist, centroids, hclusters);

					// iterate through until the desired number of clusters is reached (K)
					for (int i=0; i < (pointlist.size() - K); i++)
						maxMerge(pointlist, centroids, hclusters, dmax);

					// calculate the centroids of each group so that the SSE calculation is correct
					updateCentroids(pointlist, centroids);

					// sort lists by the group number for aesthetic reasons
					pointlist.sort(compg);
					centroids.sort(compg);

					cout << "  Completed after " << (int)(pointlist.size()-K) << " merge steps." << endl;
					cout << setprecision(3) << fixed;
					cout << "  SSE = " << SSEcalc(pointlist, centroids) << endl;
					mode = 2;

				} // end if (flag)

				break;
			case 10:  // Assemble agglomerative heirarchichal clusters using GAVG

				// checking to see if a file has already been read in, if not it reads in the selected file
				if (pointlist.size() == 0) flag = ReadInput(inputFile, pointlist, false); else flag = true;

				if (flag) {
					// assembling the initial singleton groups (each point is its own group)
					dmax = initCentroids(pointlist, centroids, hclusters);

					// iterate through until the desired number of clusters is reached (K)
					for (int i=0; i < (pointlist.size() - K); i++)
						gavgMerge(pointlist, centroids, hclusters, dmax);

					// calculate the centroids of each group so that the SSE calculation is correct
					updateCentroids(pointlist, centroids);

					// sort lists by the group number for aesthetic reasons
					pointlist.sort(compg);
					centroids.sort(compg);

					cout << "  Completed after " << (int)(pointlist.size()-K) << " merge steps." << endl;
					cout << setprecision(3) << fixed;
					cout << "  SSE = " << SSEcalc(pointlist, centroids) << endl;
					mode = 3;

				} // end if (flag)

				break;
			case 11:  // Assemble agglomerative heirarchichal clusters using CENTR

				// checking to see if a file has already been read in, if not it reads in the selected file
				if (pointlist.size() == 0) flag = ReadInput(inputFile, pointlist, false); else flag = true;

				if (flag) {
					// assembling the initial singleton groups (each point is its own group)
					dmax = initCentroids(pointlist, centroids, hclusters);

					// iterate through until the desired number of clusters is reached (K)
					for (int i=0; i < (pointlist.size() - K); i++)
						centrMerge(pointlist, centroids, hclusters, dmax);

					// calculate the centroids of each group so that the SSE calculation is correct
					updateCentroids(pointlist, centroids);

					// sort lists by the group number for aesthetic reasons
					pointlist.sort(compg);
					centroids.sort(compg);

					cout << "  Completed after " << (int)(pointlist.size()-K) << " merge steps." << endl;
					cout << setprecision(3) << fixed;
					cout << "  SSE = " << SSEcalc(pointlist, centroids) << endl;
					mode = 4;

				} // end if (flag)

				break;
			case 12:  // Display the data points by sending them to stdout
				cout << "Data Points" << endl;
				cout << "===========" << endl;
				it = pointlist.begin();
				while (it != pointlist.end()) {
					cout << *it << endl;
					it++;
				}
				break;
			case 13:   // Display the centroids by sending them to stdout
				cout << "Centroids" << endl;
				cout << "=========" << endl;
				it = centroids.begin();
				while (it != centroids.end()) {
					cout << *it << endl;
					it++;
				}
				cout << endl << "SSE = " << SSEcalc(pointlist, centroids) << endl << endl;				
				break;
			case 14:  // write data to file
				if (pointlist.size() > 0) {
					xmax = maxX(pointlist);
					xmin = minX(pointlist);
					ymax = maxY(pointlist);
					ymin = minY(pointlist);
					outdata(pointlist, xmin, xmax, ymin, ymax, K, SSEcalc(pointlist, centroids), mode, "setdata.js");
				} // end if (pointlist.size())
				else cerr << "Nothing to write." << endl;
				break;
		}

		mainmenu.draw(0,50);
		mainmenu.addenda("K = ",K,2,false);
		mainmenu.addenda("Max iterations = ",iter,5,false);
		mainmenu.addenda("Current File: " + inputFile,true);

	 } while ((idx=mainmenu.prompt()) != 0); {

	}


}