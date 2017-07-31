#ifndef POINT_H
#define POINT_H

class Point {
    private:
    //double X, Y;
    public:
    double X,Y;
    Point();
    Point(double, double);
     
    // Functions
    void read();
    void write();
    void setvalues(double&, double&);
    void getvalues(double&, double&);
    double getXvalue() const;
    double getYvalue() const;
      
};
#endif
