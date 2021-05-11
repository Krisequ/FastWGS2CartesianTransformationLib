//
// Created by Krzysztof Gromada on 10.05.2021.
//
/// Tests belowe are based on open maps readout (with relatively high inaccuracies).
/// Used mainly for code development and major error checking, more accurate comparisons should be available in the docs

#include <gtest/gtest.h>
#include "geographicLibrary.h"

point2d wptLodz     (19.426100686098042, 51.76221640076861);///< Example coordinate in Lodz city, Poland, used for testing
point2d wptBetween  (17.964919081582230, 51.73500928978871);///< Example coordinate in Poland, used for testing
point2d wptOpole    (17.992384901215985, 50.82101854251398);///< Example coordinate near Opole city, Poland, used for testing


///getAzimuthRads function tests
TEST(geographyCalculatorTests, myAzimuthEstimation){
    auto angle = rad2deg(fastGeoCalc::getAzimuthRads(wptLodz, wptOpole));
    std::cout << "Lodzi to Opole Azimuth " << angle << std::endl;
    EXPECT_NEAR(angle, 360 - 135, 2 );
    EXPECT_NEAR(rad2deg(fastGeoCalc::getAzimuthRads(wptLodz, wptBetween)), 270, 2 );
    EXPECT_NEAR(rad2deg(fastGeoCalc::getAzimuthRads(wptBetween, wptOpole)), 180, 2 );

    EXPECT_NEAR(rad2deg(fastGeoCalc::getAzimuthRads(wptOpole, wptLodz)), 45, 2 );
    EXPECT_NEAR(fabs(rad2deg(fastGeoCalc::getAzimuthRads(wptOpole, wptLodz)) - angle), 180, 2 );
}

///getRapidAzimuthRads function tests
TEST(geographyCalculatorTests, myApproxAzimuthEstimation){
    auto angle = rad2deg(fastGeoCalc::getRapidAzimuthRads(wptLodz, wptOpole));
    std::cout << "Lodzi to Opole Rapid Calc Azimuth is: " << angle << std::endl;
    EXPECT_NEAR(angle, 360 - 135, 2 );
    EXPECT_NEAR(rad2deg(fastGeoCalc::getRapidAzimuthRads(wptLodz, wptBetween)), 270, 2 );
    EXPECT_NEAR(rad2deg(fastGeoCalc::getRapidAzimuthRads(wptBetween, wptOpole)), 180, 2 );

    EXPECT_NEAR(rad2deg(fastGeoCalc::getRapidAzimuthRads(wptOpole, wptLodz)), 45, 2 );
    EXPECT_NEAR(fabs(rad2deg(fastGeoCalc::getRapidAzimuthRads(wptOpole, wptLodz)) - angle), 180, 2 );
}

///Tests of transformation from XY to UV (global coordinate system)
TEST(geographyCalculatorTests, myXYtoNECoordinates){
    double dist= 20'000, angle=135;

    /// Distance and azimuth calculation is tested above, to test XY->WGS84 the radius and angle is transformed to wgs84 and back to confirm accuracy
    ASSERT_NEAR(fastGeoCalc::getDistanceMeters(fastGeoCalc::findPointInDistanceUnderAzimuth(wptLodz, dist, angle), wptLodz),
                dist,dist*0.001);
    ASSERT_NEAR(rad2deg(fastGeoCalc::getAzimuthRads(wptLodz,
                                                    fastGeoCalc::findPointInDistanceUnderAzimuth(wptLodz, dist, angle))), angle, 1);

    dist= 200'000, angle=35;
    ASSERT_NEAR(fastGeoCalc::getDistanceMeters(fastGeoCalc::findPointInDistanceUnderAzimuth(wptLodz, dist, angle), wptLodz),
                dist,dist*0.01);
    ASSERT_NEAR(rad2deg(fastGeoCalc::getAzimuthRads(wptLodz,
                                                    fastGeoCalc::findPointInDistanceUnderAzimuth(wptLodz, dist, angle))), angle, 1.2);

    dist= 100'000, angle=360-45;
    ASSERT_NEAR(fastGeoCalc::getDistanceMeters(fastGeoCalc::findPointInDistanceUnderAzimuth(wptLodz, dist, angle), wptLodz),
                dist,dist*0.005);
    ASSERT_NEAR(rad2deg(fastGeoCalc::getAzimuthRads(wptLodz,
                                                    fastGeoCalc::findPointInDistanceUnderAzimuth(wptLodz, dist, angle))), angle, 1);

    dist= 30'000, angle=176;
    ASSERT_NEAR(fastGeoCalc::getDistanceMeters(fastGeoCalc::findPointInDistanceUnderAzimuth(wptLodz, dist, angle), wptLodz),
                dist,dist*0.001);
    ASSERT_NEAR(rad2deg(fastGeoCalc::getAzimuthRads(wptLodz,
                                                    fastGeoCalc::findPointInDistanceUnderAzimuth(wptLodz, dist, angle))), angle, 1);

}







