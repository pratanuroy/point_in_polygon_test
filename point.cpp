#include <iostream>
#include "point.h"

Point::Point() {
X = 0.0;
Y = 0.0;
}
Point::Point(double posX, double posY) {
        X = posX;
        Y = posY;
}
   
void Point::write() {
std::cout << "<Point(" << X << "," << Y << ")>" << std::endl;
}

void Point::read(){
std::cin >> X;
std::cin >> Y;
}

void Point::setvalues(double &posX, double &posY)
{
    X = posX;
    Y = posY;
}

void Point::getvalues(double &posX, double &posY)
{
    posX = X;
    posY = Y;
}

double Point::getXvalue() const
{
    return X;
}
double Point::getYvalue() const
{
    return Y;
}
