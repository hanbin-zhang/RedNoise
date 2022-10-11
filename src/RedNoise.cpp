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
#include <TexturePoint.h>
#include <TextureMap.h>
#include "glm/mat3x3.hpp"

#define WIDTH 320
#define HEIGHT 240

struct {
    CanvasPoint from;
    CanvasPoint to;
} line;

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

void draw_texture_line(DrawingWindow &window, CanvasPoint from, CanvasPoint to, TextureMap textureMap) {
    float x_diff = to.x - from.x;
    float y_diff = to.y - from.y;

    float numberOfSteps = std::max(abs(x_diff), abs(y_diff));
    float x_step_size = x_diff / numberOfSteps;
    float y_step_size = y_diff / numberOfSteps;

    float texture_x_diff = to.texturePoint.x - from.texturePoint.x;
    float texture_y_diff = to.texturePoint.y - from.texturePoint.y;

    //float numberOfStepsTexture = std::max(abs(texture_x_diff), abs(texture_y_diff));
    float x_step_size_Texture = texture_x_diff / numberOfSteps;
    float y_step_size_Texture = texture_y_diff / numberOfSteps;

    for (float i = 0.0; i < numberOfSteps; ++i) {
        float x = from.x + i*x_step_size;
        float y = from.y + i*y_step_size;

        float texture_x = from.texturePoint.x + i*x_step_size_Texture;
        float texture_y = from.texturePoint.y + i*y_step_size_Texture;

        int texture_x_int = int (texture_x) % int(textureMap.width);
        int texture_y_int = int (texture_y) % int(textureMap.height);

        window.setPixelColour(size_t(round(x)), size_t(round(y)),
                              textureMap.pixels[textureMap.width*texture_y_int-(textureMap.width-texture_x_int)]);
    }
}

/*std::array<float, 2> calculate_step_size(std::array<float, 2>) {
    float x_diff = to.x - from.x;
    float y_diff = to.y - from.y;

    float numberOfSteps = std::max(abs(x_diff), abs(y_diff));
    float x_step_size = x_diff / numberOfSteps;
    float y_step_size = y_diff / numberOfSteps;

    return std::array<float, 2> {x_step_size, y_step_size};
}*/

void fill_half_triangle(DrawingWindow &window, CanvasPoint from_start, CanvasPoint to_start, CanvasPoint from_end, CanvasPoint to_end, const Colour& colour) {

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

std::vector<CanvasPoint> interpolatingCanvasPoint(CanvasPoint from, CanvasPoint to, float step_number) {
    float x_diff = from.x - to.x;
    float y_diff = from.y - to.y;

    float x_step_size = x_diff / step_number;
    float y_step_size = y_diff/ step_number;

    float x_diff_texture = from.texturePoint.x - to.texturePoint.x;
    float y_diff_texture = from.y - to.y;

    float x_step_size_texture = x_diff_texture / step_number;
    float y_step_size_texture = y_diff_texture / step_number;

    std::vector<CanvasPoint> result;
    for (float i = 0.0; i < step_number; ++i) {
        float x_from = from.x + i * x_step_size;
        float y_from = from.y + i * y_step_size;

        float x_texture = from.texturePoint.x + i * x_step_size_texture;
        float y_texture = to.texturePoint.y + i * y_step_size_texture;

        CanvasPoint canvasPoint = CanvasPoint(x_from, y_from);

        from.texturePoint = TexturePoint(x_texture, y_texture);
        to.texturePoint = TexturePoint(x_texture, y_texture);

        result.push_back(canvasPoint);
    }
    return result;
}

void fill_half_texture_triangle(DrawingWindow &window, CanvasPoint from_start, CanvasPoint to_start, CanvasPoint from_end, CanvasPoint to_end, TextureMap textureMap) {

    float x_diff_to = from_end.x - to_start.x;
    float y_diff_to = from_end.y - to_start.y;

    float x_diff_from = to_end.x - from_start.x;
    float y_diff_from = to_end.y - from_start.y;

    float x_diff_to_texture = from_end.texturePoint.x - to_start.texturePoint.x;
    float y_diff_to_texture = from_end.y - to_start.y;

    float x_diff_from_texture = to_end.texturePoint.x - from_start.texturePoint.x;
    float y_diff_from_texture = to_end.texturePoint.y - from_start.texturePoint.y;

    float numberOfSteps_x = std::max(abs(x_diff_from), abs(x_diff_to));
    float numberOfSteps_y = std::max(abs(y_diff_from), abs(y_diff_to));
    float numberOfSteps_x_texture = std::max(abs(x_diff_from_texture), abs(x_diff_to_texture));
    float numberOfSteps_y_texture = std::max(abs(y_diff_from_texture), abs(y_diff_to_texture));
    float numberOfSteps1 = std::max(numberOfSteps_x, numberOfSteps_y);
    float numberOfSteps2 = std::max(numberOfSteps_x_texture, numberOfSteps_y_texture);
    float numberOfSteps = std::max(numberOfSteps1, numberOfSteps2);

    float x_step_size_from = x_diff_from / numberOfSteps;
    float y_step_size_from = y_diff_from / numberOfSteps;
    float x_step_size_to = x_diff_to / numberOfSteps;
    float y_step_size_to = y_diff_to / numberOfSteps;



    float x_step_size_from_texture = x_diff_from_texture / numberOfSteps;
    float y_step_size_from_texture = y_diff_from_texture / numberOfSteps;
    float x_step_size_to_texture = x_diff_to_texture / numberOfSteps;
    float y_step_size_to_texture = y_diff_to_texture / numberOfSteps;

    for (float i = 0.0; i <= numberOfSteps; ++i) {
        float x_from = from_start.x + i * x_step_size_from;
        float y_from = from_start.y + i * y_step_size_from;

        float x_to = to_start.x + i * x_step_size_to;
        float y_to = to_start.y + i * y_step_size_to;

        float x_from_texture = from_start.texturePoint.x + i * x_step_size_from_texture;
        float y_from_texture = from_start.texturePoint.y + i * y_step_size_from_texture;

        float x_to_texture = to_start.texturePoint.x + i * x_step_size_to_texture;
        float y_to_texture = to_start.texturePoint.y + i * y_step_size_to_texture;

        CanvasPoint from = CanvasPoint(x_from, y_from);
        CanvasPoint to = CanvasPoint(x_to, y_to);

        from.texturePoint = TexturePoint(x_from_texture, y_from_texture);
        to.texturePoint = TexturePoint(x_to_texture, y_to_texture);

        draw_texture_line(window,from, to, textureMap);
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
    fill_half_triangle(window,  vertices[2], vertices[2], mid_point, vertices[1],colour);
}

void draw_stroked_triangles(DrawingWindow &window, CanvasTriangle triangle, const Colour& colour) {
    draw_line(window, triangle.v0(), triangle.v1(), colour);
    draw_line(window, triangle.v1(), triangle.v2(), colour);
    draw_line(window, triangle.v2(), triangle.v0(), colour);
}

void draw_filled_triangles(DrawingWindow &window, CanvasTriangle triangle, const Colour& colour) {
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

void textureMapper(DrawingWindow &window, CanvasTriangle canvasTriangle, TextureMap textureMap) {
    std::array<CanvasPoint, 3> vertices = canvasTriangle.vertices;
    while (true) {
        if (vertices[0].y <= vertices[1].y && vertices[1].y <= vertices[2].y) break;
        else if (vertices[2].y <= vertices[0].y) std::swap(vertices[2], vertices[0]);
        else if (vertices[1].y <= vertices[0].y) std::swap(vertices[0], vertices[1]);
        else if (vertices[2].y <= vertices[1].y) std::swap(vertices[2], vertices[1]);
    }

    draw_stroked_triangles(window, canvasTriangle, Colour(255, 255, 255));
    CanvasPoint midpoint = find_mid_point(vertices);
    std::array<CanvasPoint, 3> texture_vertices{
        CanvasPoint(vertices[0].texturePoint.x, vertices[0].texturePoint.y),
        CanvasPoint(vertices[1].texturePoint.x, vertices[1].texturePoint.y),
        CanvasPoint(vertices[2].texturePoint.x, vertices[2].texturePoint.y)
    };
    CanvasPoint mid_texture_point = find_mid_point(texture_vertices);
    std::cout << vertices[0] << std::endl;
    std::cout << vertices[1] << std::endl;
    std::cout << vertices[2] << std::endl;

    midpoint.texturePoint = TexturePoint(mid_texture_point.x, mid_texture_point.y);

    fill_half_texture_triangle(window, vertices[0], vertices[0], midpoint, vertices[1], textureMap);
    fill_half_texture_triangle(window,  vertices[2], vertices[2],midpoint, vertices[1], textureMap);

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

glm::mat3x3 reverse_mtx(glm::mat3x3 mat) {
    glm::mat3x3 reverse_matrix;
    int determinant;
    for(int i = 0; i < 3; i++) {
        determinant = determinant + (mat[0][i] * (mat[1][(i + 1) % 3] * mat[2][(i + 2) % 3] -
                                                  mat[1][(i + 2) % 3] * mat[2][(i + 1) % 3]));
        std::cout << determinant << std::endl;
    }
    std::cout << determinant << std::endl;
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            reverse_matrix[i][j] = ((mat[(j+1)%3][(i+1)%3] * mat[(j+2)%3][(i+2)%3]) - (mat[(j+1)%3][(i+2)%3] * mat[(j+2)%3][(i+1)%3])) / determinant;
        }
    }
    return reverse_matrix;
}

glm::mat3x3 calculate_affine_mtx(CanvasTriangle triangle) {
    glm::mat3x3 sdl_matrix = glm::mat3x3(triangle.v0().x, triangle.v1().x, triangle.v2().x,
                                         triangle.v0().y, triangle.v1().y, triangle.v2().y,
                                         1.0, 1.0, 1.0);
    glm::mat3x3 texture_matrix = glm::mat3x3(triangle.v0().texturePoint.x, triangle.v1().texturePoint.x, triangle.v2().texturePoint.x,
                                         triangle.v0().texturePoint.y, triangle.v1().texturePoint.y, triangle.v2().texturePoint.y,
                                         1.0, 1.0, 1.0);
    glm::mat3x3 sdl_mtx_inverse;



    return reverse_mtx(sdl_matrix);
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
    std::vector<CanvasTriangle> triangles;
    std::vector<Colour> colours;
    //add_triangle(triangles, colours);

    CanvasPoint v0 = CanvasPoint(160, 10);
    CanvasPoint v1 = CanvasPoint(300, 230);
    CanvasPoint v2 = CanvasPoint(10, 150);
    v0.texturePoint = TexturePoint(195, 5);
    v1.texturePoint = TexturePoint(395, 380);
    v2.texturePoint = TexturePoint(65, 330);
    CanvasTriangle texture_triangle = CanvasTriangle(v0, v1, v2);

    TextureMap textureMap = TextureMap("texture.ppm");
    //textureMapper(window,  texture_triangle, textureMap);
    glm::mat3x3 matrix = calculate_affine_mtx(texture_triangle);
    std::cout << matrix[0][0] << ", " << matrix[0][1] << ", "  << matrix[0][2] << std::endl;
    std::cout << matrix[1][0] << ", "<< matrix[1][1] << ", "  << matrix[1][2] << std::endl;
    std::cout << matrix[2][0]<< ", " << matrix[2][1]   << ", "<< matrix[2][2] << std::endl;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window, triangles, colours);
		//draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !

        /*for(int i = 0; i < int(triangles.size()); ++i) {
            draw_filled_triangles(window, triangles[i], colours[i]);
        }*/

        window.renderFrame();
        //break;
	}
}
