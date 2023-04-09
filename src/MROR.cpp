#include "MROR.h"

void MROR::SetRadius(double radius_multiplier) {
    radius_multiplier_ = radius_multiplier;
}

double MROR::GetRadius() {
    return radius_multiplier_;
}

void MROR::SetAngle(double azimuth_angle) {
    azimuth_angle_ = azimuth_angle;
}

double MROR::GetAngle() {
    return azimuth_angle_;
}

void MROR::SetMinNeighbors(double min_neighbors) {
    min_neighbors_ = min_neighbors;
}

double MROR::GetMinNeighbors() {
    return min_neighbors_;
}

void MROR::SetMinRadius(double min_search_radius) {
    min_search_radius_ = min_search_radius;
}

double MROR::GetMinSearchRadius() {
    return min_search_radius_;
}

void MROR::SetKSearch(double KSearch){
    KSearch_ = KSearch;
}

double MROR::GetKSearch(){
    return KSearch_;
}
