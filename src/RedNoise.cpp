#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include "glm/vec3.hpp"
#include <CanvasPoint.h>
#include <Colour.h>
#include <CanvasTriangle.h>
#include <stdlib.h>
#include <map>

#define WIDTH 320
#define HEIGHT 240

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {

    std::vector<float> result;
    float  step = (to - from) / float(numberOfValues-1);
    for (int i = 0; i < numberOfValues; i++) {
        result.push_back(from + float(i)*step);
    }
    return result;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues) {
    std::vector<glm::vec3> resultV3;

    std::vector<float> result1;
    std::vector<float> result2;
    std::vector<float> result3;

    result1 = interpolateSingleFloats(from.x, to.x, numberOfValues);
    result2 = interpolateSingleFloats(from.y, to.y, numberOfValues);
    result3 = interpolateSingleFloats(from.z, to.z, numberOfValues);

    for (int i = 0; i < numberOfValues; i++) {
        resultV3.emplace_back(result1[i], result2[i], result3[i]);
    }
    return resultV3;
}

uint32_t colour_uint32(const Colour& colour) {
    return (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
}

void draw_line(DrawingWindow &window, CanvasPoint from, CanvasPoint to, const Colour& colour) {
    float x_diff = to.x - from.x;
    float y_diff = to.y - from.y;

    float numberOfSteps = std::max(abs(x_diff), abs(y_diff));
    float x_step_size = x_diff / numberOfSteps;
    float y_step_size = y_diff / numberOfSteps;

    for (float i = 0.0; i < numberOfSteps; ++i) {
        float x = from.x + i*x_step_size;
        float y = from.y + i*y_step_size;

        window.setPixelColour(size_t(round(x)), size_t(round(y)), colour_uint32(colour));
    }
}

void fill_half_triangle(DrawingWindow &window, CanvasPoint from_start, CanvasPoint to_start, CanvasPoint from_end, CanvasPoint to_end, const Colour& colour) {
    /*CanvasPoint from_start = vertices[0];
    CanvasPoint to_start = vertices[0];

    CanvasPoint from_end = mid_point;
    CanvasPoint to_end = vertices[1];*/

    float x_diff_to = from_end.x - to_start.x;
    float y_diff_to = from_end.y - to_start.y;

    float x_diff_from = to_end.x - from_start.x;
    float y_diff_from = to_end.y - from_start.y;

    float numberOfSteps_x = std::max(abs(x_diff_from), abs(x_diff_to));
    float numberOfSteps_y = std::max(abs(y_diff_from), abs(y_diff_to));
    float numberOfSteps = std::max(numberOfSteps_x, numberOfSteps_y);

    float x_step_size_from = x_diff_from / numberOfSteps;
    float y_step_size_from = y_diff_from / numberOfSteps;
    float x_step_size_to = x_diff_to / numberOfSteps;
    float y_step_size_to = y_diff_to / numberOfSteps;

    for (float i = 0.0; i < numberOfSteps; ++i) {
        float x_from = from_start.x + i * x_step_size_from;
        float y_from = from_start.y + i * y_step_size_from;

        float x_to = to_start.x + i * x_step_size_to;
        float y_to = to_start.y + i * y_step_size_to;

        draw_line(window, CanvasPoint(x_from, y_from), CanvasPoint(x_to, y_to), colour);
    }
}

CanvasPoint find_mid_point(std::array<CanvasPoint, 3> vertices) {
    CanvasPoint mid_point;
    mid_point.y = vertices[1].y;

    auto top_bottom_x_diff = float(abs(int(vertices[0].x - vertices[2].x)));
    auto top_bottom_y_diff = float (abs(int(vertices[0].y - vertices[2].y)));

    float mid_x_diff = top_bottom_x_diff/top_bottom_y_diff * float(abs(int(vertices[0].y - vertices[1].y)));

    if (vertices[0].x >= vertices[2].x) mid_point.x = vertices[0].x - mid_x_diff;
    else mid_point.x = vertices[0].x + mid_x_diff;
    return mid_point;
}

void fill_triangle(DrawingWindow &window, CanvasTriangle triangle, const Colour& colour) {
    std::array<CanvasPoint, 3> vertices = triangle.vertices;
    while (true) {
        if (vertices[0].y <= vertices[1].y && vertices[1].y <= vertices[2].y) break;
        else if (vertices[2].y <= vertices[0].y) std::swap(vertices[2], vertices[0]);
        else if (vertices[1].y <= vertices[0].y) std::swap(vertices[0], vertices[1]);
        else if (vertices[2].y <= vertices[1].y) std::swap(vertices[2], vertices[1]);
    }

    CanvasPoint mid_point = find_mid_point(vertices);

    // draw top triangle
    fill_half_triangle(window, vertices[0], vertices[0], mid_point, vertices[1], colour);
    // draw bottom triangle
    fill_half_triangle(window, mid_point, vertices[1], vertices[2], vertices[2], colour);
}

void draw_stroked_triangles(DrawingWindow &window, CanvasTriangle triangle, const Colour& colour) {
    draw_line(window, triangle.v0(), triangle.v1(), colour);
    draw_line(window, triangle.v1(), triangle.v2(), colour);
    draw_line(window, triangle.v2(), triangle.v0(), colour);
    fill_triangle(window, triangle, colour);
}

void draw(DrawingWindow &window) {
	window.clearPixels();

    glm::vec3 topLeft(255, 0, 0);        // red
    glm::vec3 topRight(0, 0, 255);       // blue
    glm::vec3 bottomRight(0, 255, 0);    // green
    glm::vec3 bottomLeft(255, 255, 0);   // yellow

    std::vector<glm::vec3> left = interpolateThreeElementValues(topLeft, bottomLeft, int(window.width));
    std::vector<glm::vec3> right = interpolateThreeElementValues(topRight, bottomRight, int(window.width));

	for (size_t y = 0; y < window.height; y++) {
        std::vector<glm::vec3> current_line = interpolateThreeElementValues(left[int (y)], right[int (y)], int(window.width));
		for (size_t x = 0; x < window.width; x++) {
			uint32_t colour = (255 << 24) + (int(current_line[int (x)].x) << 16)
                    + (int(current_line[int (x)].y) << 8) + int(current_line[int (x)].z);
			window.setPixelColour(x, y, colour);
		}
	}
}

void add_triangle(std::vector<CanvasTriangle>& triangles, std::vector<Colour>& colours) {
    Colour colour = Colour(rand()%255+1, rand()%255+1, rand()%255+1);
    CanvasTriangle triangle = CanvasTriangle(
            CanvasPoint(float(rand()%WIDTH+1), float(rand()%HEIGHT+1)),
            CanvasPoint(float(rand()%WIDTH+1), float(rand()%HEIGHT+1)),
            CanvasPoint(float(rand()%WIDTH+1), float(rand()%HEIGHT+1)));
    triangles.push_back(triangle);
    colours.push_back(colour);
}

void handleEvent(SDL_Event event, DrawingWindow &window, std::vector<CanvasTriangle>& triangles, std::vector<Colour>& colours) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
        else if (event.key.keysym.sym == SDLK_u) add_triangle(triangles, colours);
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
    std::vector<CanvasTriangle> triangles;
    std::vector<Colour> colours;
    //add_triangle(triangles, colours);
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window, triangles, colours);
		//draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !

        for(int i = 0; i < int(triangles.size()); ++i) {
            draw_stroked_triangles(window, triangles[i], colours[i]);
        }

        window.renderFrame();
        //break;
	}
}
