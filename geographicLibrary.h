//
// Created by Krzysztof Gromada on 10.05.2021.
//

#ifndef GEOGRAPHICLIBRARY_H
#define GEOGRAPHICLIBRARY_H
#include <cmath>

constexpr double EarthRadius=6'367'450; ///< Average earth radius (according to WGS84)- in meters 6378100+6356800 /2
constexpr double MetersPerLatitudeDegree=6'367'450*2*M_PI/360; ///< Meters per degree approximation
constexpr double LatitudeDegreePerMeter=1/MetersPerLatitudeDegree; ///< Degrees per meter approximation

inline double rad2deg(const double &rad) {
    constexpr double mult = 180.0 / M_PI;
    return rad * mult;
}
inline double deg2rad(const double &deg) {
    constexpr double mult =  M_PI/180.0;
    return deg * mult;
}

struct point2d{
    double x, y;
    point2d(const point2d &in):x(in.x),y(in.y){}
    explicit point2d(const double &x=NAN, const double &y=NAN):x(x),y(y){}
};



namespace fastGeoCalc{
    /// Returns azimuth angle from ptIn to ptOut. Assumption of spheroid earth is taken
    ///
    /// eauation taken   tan(azimuth) = sin(L) / (cos(p1)&tan(p2)-sin(p1)cos(L),
    /// where ptIn=(Lon = 0,Lat = p1), ptOut = (Lon = L, Lat = p2)
    ///
    /// so the following equation can be deduced:
    ///
    ///     azimuth = atan( sin(l2-l1) / (cos(p1)&tan(p2)-sin(p1)cos(l2-l1) )
    ///
    /// \param ptIn - angular position on earth (degrees)
    /// \param ptOut - angular position on earth (degrees)
    /// \return float angle from north (RADIANS)
    inline double getAzimuthRads(const point2d &ptIn, const point2d &ptOut){
        double l1 = deg2rad(ptIn.x),
        p1 = deg2rad( ptIn.y  ),
        l2 = deg2rad( ptOut.x ),
        p2 = deg2rad( ptOut.y );

        double output = atan2(sin(l2-l1), cos(p1)*tan(p2)-sin(p1)*cos(l2-l1));

        if(output<0){/// <-pi, pi> to <0, 2pi>
            return 2*M_PI+output;
        }

        return output;
    }
    /// Return approximate azimuth angle between 2 points on earth surface - accurate for short distances
    /// Using this algorithm: http://zasoby.open.agh.edu.pl/~11sjjurek/azymut.html
    /// \param ptIn - angular position on earth (degrees)
    /// \param ptOut - angular position on earth (degrees)
    /// \return double angle from north (RADIANS)  - can be change to float due to low accuracy
    inline double getRapidAzimuthRads(const point2d &ptIn, const point2d &ptOut){
        auto dy = (ptOut.y-ptIn.y);
        auto dx = ((ptOut.x-ptIn.x)*cos(deg2rad((ptIn.y+ptOut.y)/2))); //scaling the X axis
        double output = atan2(dx,dy); // not dy, dx due to different coordinate system (angle from y not x)
        if(output<0){
            return 2*M_PI+output; /// <-pi, pi> to <0, 2pi>
        }

        return output;
    }
    /// Interface for Haversine formula - returns distance between 2 points on the ground
    /// \param ptIn - angular position on earth (degrees)
    /// \param ptOut - angular position on earth (degrees)
    /// \return distance in meters
    inline double getDistanceMeters(const point2d &ptIn, const point2d &ptOut){
        double lat1 = deg2rad(ptIn.y), lat2=deg2rad(ptOut.y);
        return 2.*EarthRadius*(sqrt(
                pow(sin((lat2-lat1)/2.),2) + cos(lat2)*cos(lat1)*pow(sin(deg2rad(ptOut.x-ptIn.x)/2.),2)
        )); //używany jest zapis oparty o sinusy a nie cosinusy w celu zmniejszenie błędów numerycznych
    }
    // /////////////////////// Geo to XYZ ///////////////////////

    /// Funkcja pozwalająca na zmianę układu współrzędnych punktu.
    /// pozycji NW na kołowy układ współrzędnych XY (opis w postaci kątu azymutu i odległości). Pierwszy punkt określa środek układu współrzędnych.
    ///
    /// Zgodnie z testami dla 100km odległości dokładność powinna być na poziomie +/- 500m. Dla 200km może przekraczać 1.5km.
    ///
    /// @param zeroWptDegs - środek układu współrzędnych względem którego będą liczone pozycje (najpewniej wsp. samolotu)
    /// @param wptDegs - współrzędne punktu rzucanego
    /// @return zwraca punkt we współrzędnych XY
    inline point2d movePointToXYcs(const point2d &zeroWptDegs, const point2d &wptDegs){
        double azimuthAngleDegs = getAzimuthRads(zeroWptDegs, wptDegs);
        double r = getDistanceMeters(zeroWptDegs,wptDegs);
        return point2d(sin(rad2deg(azimuthAngleDegs))*r,cos(rad2deg(azimuthAngleDegs))*r);
    }

    // /////////////////////// XYZ to Geo ///////////////////////

    /// Reverse function to "movePointToXYcs" changes XY position to global coordinate system.
    inline point2d findPointInDistanceUnderAzimuth (const point2d& startingWptDegs, const double groundDistanceInMeters, double azimuthInDegs){
        double  dx = groundDistanceInMeters*sin(deg2rad(azimuthInDegs))*LatitudeDegreePerMeter/cos(deg2rad(startingWptDegs.y)),
                dy = groundDistanceInMeters*cos(deg2rad(azimuthInDegs))*LatitudeDegreePerMeter;
        return point2d(startingWptDegs.x+dx, startingWptDegs.y+dy);
    }
}






/////////////////////// For Lighter version comment all the code below ///////////////////////

#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/strategies/spherical/distance_haversine.hpp>
#include <boost/geometry/algorithms/distance.hpp>

namespace  bg = boost::geometry;
typedef  bg::cs::geographic<bg::degree> geoCS;
typedef  bg::strategy::distance::haversine<double> stratLinDist;
typedef  bg::strategy::azimuth::geographic<> stratAzimuth;


namespace boostInterface {
    /// Boost::Geometry equivalent to the function proposed above
    inline double getAzimuthRads(const bg::model::point<double,2,geoCS> &ptIn, const bg::model::point<double,2,geoCS> &ptOut){
        double output = 0, output2 = 0;
        stratAzimuth().apply( deg2rad(bg::get<0>(ptIn)),
                            deg2rad(bg::get<1>(ptIn)),
                            deg2rad(bg::get<0>(ptOut)),
                            deg2rad(bg::get<1>(ptOut)),
                            output,output2);
        return output;
    }

    /// Boost::Geometry equivalent to the function proposed above
    inline double getDistanceMeters(const bg::model::point<double,2,geoCS> &ptIn, const bg::model::point<double,2,geoCS> &ptOut){
        return bg::distance(ptIn, ptOut, stratLinDist()) * EarthRadius ;
    }
};




#endif //GEOGRAPHICLIBRARY_H
