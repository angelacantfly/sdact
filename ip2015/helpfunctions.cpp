
//
//  helpfunctions.cpp
//  ip2015
//
//  Created by ANGELA ZHOU on 5/3/15.
//  Copyright (c) 2015 Z Sweedyk. All rights reserved.
//

#include "helpfunctions.h"
vector<Point> find_boundary(int width, int height){
    vector<Point> boundaryPoints;
    int center_x = width/2;
    int center_y = height/2;
    int square_size = (center_x < center_y? center_x:center_y);
    
    Point top_left, top_right, bot_left, bot_right;
    top_left.x = center_x  - (square_size/2);
    top_left.y = center_y - (square_size/2);
    top_right.x = top_left.x + square_size;
    top_right.y = top_left.y;
    
    bot_left.x = top_left.x;
    bot_left.y = top_left.y + square_size;
    bot_right.x = top_right.x;
    bot_right.y = bot_left.y;
    
    for (int i = top_right.x; i > top_left.x; i--)
    {
        Point newPoint;
        newPoint.x =i;
        newPoint.y = top_left.y;
        boundaryPoints.push_back(newPoint);
    }
    
    for (int i = top_left.y; i < bot_left.y ; i++) {
        Point newPoint;
        newPoint.x = top_left.x;
        newPoint.y = i;
        boundaryPoints.push_back(newPoint);
    }
    
    for (int i = bot_left.x; i < bot_right.x; i++) {
        Point newPoint;
        newPoint.x = i;
        newPoint.y = bot_left.y;
        boundaryPoints.push_back(newPoint);
    }
    for (int i = bot_right.y; i >top_right.y; i--) {
        Point newPoint;
        newPoint.x = top_right.x;
        newPoint.y = i;
        boundaryPoints.push_back(newPoint);
    }
    
//    boundaryPoints.pop_back();
    cout << "boundary point size : " <<boundaryPoints.size() << endl;
    return boundaryPoints;
}

void find_vertices_square(int width, int height, Point *top_left, Point *top_right, Point *bot_left, Point *bot_right){
    int center_x = width/2;
    int center_y = height/2;
    int square_size = (center_x > center_y? center_x:center_y);
    
    top_left->x = center_x  - (square_size/2);
    top_left->y = center_y - (square_size/2);
    top_right->x = top_left->x + square_size;
    top_right->y = top_left->y;
    
    bot_left->x = top_left->x;
    bot_left->y = top_left->y + square_size;
    bot_right->x = top_right->x;
    bot_right->y = bot_left->y;
    
}