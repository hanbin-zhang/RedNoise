#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <vector>
#include <algorithm>
#include "glm/vec3.hpp"
#include <CanvasPoint.h>
#include <Colour.h>
#include <stdlib.h>
#include <map>
#include <TexturePoint.h>
#include <TextureMap.h>
#include "glm/mat3x3.hpp"
#include <iostream>
#include "ModelTriangle.h"

#define WIDTH 640
#define HEIGHT 480

float depth_buffer[WIDTH][HEIGHT];

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

std::map<std::string, Colour> read_colour_palette(const std::string& file_name) {
    std::string current_line;
    std::ifstream MyReadFile(file_name);
    std::string current_colour_string;
    std::map<std::string, Colour> colour_map;

    while (getline (MyReadFile, current_line)) {
        if (current_line.compare(0, 6, "newmtl")==0) {
            current_colour_string = split(current_line, ' ')[1];
        } else if (current_line.compare(0, 2, "Kd")==0) {
            std::vector<std::string> rgb_values = split(current_line, ' ');
            colour_map[current_colour_string] = Colour(current_colour_string,
                                                       int(std::stof(rgb_values[1])*255.0),
                                                       int(std::stof(rgb_values[2])*255.0),
                                                       int(std::stof(rgb_values[3])*255.0));
        }
    }
    return colour_map;
}

std::vector<ModelTriangle> read_OBJ_files(const std::string& file_name,
                                          const std::string& colour_file_name,
                                          float scaling) {

    std::string current_line;
    std::ifstream MyReadFile(file_name);
    std::vector<glm::vec3> vertices;
    std::vector<ModelTriangle> triangles;
    Colour current_colour;
    std::map<std::string, Colour> colour_palette = read_colour_palette(colour_file_name);

    while (getline (MyReadFile, current_line)) {
        if (current_line[0] == 'v') {
            //put a vertex on to the list
            std::vector<std::string> vertices_string = split(current_line, ' ');
            vertices.emplace_back(std::stof(vertices_string[1]) * scaling,
                                  std::stof(vertices_string[2]) * scaling,
                                  std::stof(vertices_string[3]) * scaling);

        } else if (current_line.compare(0, 6, "usemtl")==0) {
            current_colour = colour_palette[split(current_line, ' ')[1]];
        } else if (current_line[0] == 'f') {
            // add a facet triangle
            std::vector<std::string> facets_string = split(current_line, ' ');
            ModelTriangle current_triangle;
            current_triangle.colour = current_colour;
            for (int i = 1; i < int(facets_string.size()); ++i) {
                std::vector<std::string> facets_vertex_string = split(facets_string[i], '/');

                current_triangle.vertices[i-1] = vertices[std::stoi(facets_vertex_string[0])-1];
            }

            triangles.push_back(current_triangle);
        }
    }

    MyReadFile.close();
    return triangles;
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

void draw_texture_line(DrawingWindow &window, CanvasPoint from, CanvasPoint to, TextureMap textureMap, glm::mat3x3 affineMtx) {

    float texture_x_diff = to.texturePoint.x - from.texturePoint.x;
    float texture_y_diff = to.texturePoint.y - from.texturePoint.y;
    float numberOfSteps = std::max(abs(texture_x_diff), abs(texture_y_diff));
    //float numberOfStepsTexture = std::max(abs(texture_x_diff), abs(texture_y_diff));
    float x_step_size_Texture = texture_x_diff / numberOfSteps;
    float y_step_size_Texture = texture_y_diff / numberOfSteps;

    for (float i = 0.0; i < numberOfSteps; ++i) {

        float texture_x = from.texturePoint.x + i*x_step_size_Texture;
        float texture_y = from.texturePoint.y + i*y_step_size_Texture;

        int texture_x_int = int (texture_x) % int(textureMap.width);
        int texture_y_int = int (texture_y) % int(textureMap.height);

        float x = texture_x*affineMtx[0][0] + texture_y*affineMtx[0][1] + affineMtx[0][2];
        float y = texture_x*affineMtx[1][0] + texture_y*affineMtx[1][1] + affineMtx[1][2];

        window.setPixelColour(size_t(round(x)), size_t(round(y)),
                              textureMap.pixels[textureMap.width*texture_y_int-(textureMap.width-texture_x_int)]);
    }
}

void draw_line_with_depth(DrawingWindow &window, CanvasPoint from, CanvasPoint to, const Colour& colour) {
    float x_diff = to.x - from.x;
    float y_diff = to.y - from.y;

    float numberOfSteps = std::max(abs(x_diff), abs(y_diff));
    float x_step_size = x_diff / numberOfSteps;
    float y_step_size = y_diff / numberOfSteps;

    std::vector<float> depths = interpolateSingleFloats(from.depth, to.depth, int(round(numberOfSteps+1)));

    for (float i = 0; i < numberOfSteps; ++i) {
        float x = from.x + i*x_step_size;
        float y = from.y + i*y_step_size;
        if (1/depths[int (i)]>depth_buffer[int(round(x))][int(round(y))]) {
            window.setPixelColour(size_t(round(x)), size_t(round(y)), colour_uint32(colour));
            depth_buffer[int(round(x))][int(round(y))] = 1/depths[int (i)];
        }
        std::cout << 1/depths[int (i)] << ". " << depth_buffer[int(round(x))][int(round(y))] << std::endl;
        //window.setPixelColour(size_t(round(x)), size_t(round(y)), colour_uint32(colour));
    }

}

void fill_half_triangle(DrawingWindow &window, CanvasPoint start,
                        CanvasPoint from_end, CanvasPoint to_end, const Colour& colour) {
    float x_diff_to = from_end.x - start.x;
    float x_diff_from = to_end.x - start.x;
    float y_diff = to_end.y - start.y;

    float numberOfSteps_x = std::max(abs(x_diff_from), abs(x_diff_to));
    float numberOfSteps = std::max(numberOfSteps_x, abs(y_diff));

    float x_step_size_from = x_diff_from / numberOfSteps;
    float x_step_size_to = x_diff_to / numberOfSteps;
    float y_step_size = y_diff / numberOfSteps;

    std::vector<float> depth_to = interpolateSingleFloats(start.depth, to_end.depth, int(numberOfSteps));
    std::vector<float> depth_from = interpolateSingleFloats( start.depth, from_end.depth,int(numberOfSteps));


    for (float i = 0.0; i <= numberOfSteps; ++i) {
        float x_from = start.x + i * x_step_size_from;
        float y_from = start.y + i * y_step_size;

        float x_to = start.x + i * x_step_size_to;
        float y_to = start.y + i * y_step_size;
        draw_line_with_depth(window, CanvasPoint(x_from, y_from, depth_from[int(i)]),
                             CanvasPoint(x_to, y_to, depth_to[int(i)]),
                             colour);
    }
}



std::vector<CanvasPoint> interpolatingCanvasPoint(CanvasPoint from, CanvasPoint to, float stepSize) {
    float x_diff = from.x - to.x;
    float y_diff = from.y - to.y;

    float x_step_size = x_diff / stepSize;
    float y_step_size = y_diff / stepSize;

    std::vector<CanvasPoint> result;

    for (float i = 0.0; i < stepSize; ++i) {
        float x_from = from.x + i * x_step_size;
        float y_from = from.y + i * y_step_size;

        CanvasPoint canvasPoint = CanvasPoint(x_from, y_from);



        result.push_back(canvasPoint);
    }
    return result;
}

void fill_half_texture_triangle(DrawingWindow &window,
                                CanvasPoint from_start, CanvasPoint to_start,
                                CanvasPoint from_end, CanvasPoint to_end,
                                TextureMap textureMap, glm::mat3x3 affine_matrix) {

    float x_diff_to_texture = from_end.texturePoint.x - to_start.texturePoint.x;
    float y_diff_to_texture = from_end.texturePoint.y - to_start.texturePoint.y;

    float x_diff_from_texture = to_end.texturePoint.x - from_start.texturePoint.x;
    float y_diff_from_texture = to_end.texturePoint.y - from_start.texturePoint.y;

    float numberOfSteps_x_texture = std::max(abs(x_diff_from_texture), abs(x_diff_to_texture));
    float numberOfSteps_y_texture = std::max(abs(y_diff_from_texture), abs(y_diff_to_texture));
    float numberOfSteps = std::max(numberOfSteps_x_texture, numberOfSteps_y_texture);

    float x_step_size_from_texture = x_diff_from_texture / numberOfSteps;
    float y_step_size_from_texture = y_diff_from_texture / numberOfSteps;
    float x_step_size_to_texture = x_diff_to_texture / numberOfSteps;
    float y_step_size_to_texture = y_diff_to_texture / numberOfSteps;

    for (float i = 0.0; i <= numberOfSteps; ++i) {

        float x_from_texture = from_start.texturePoint.x + i * x_step_size_from_texture;
        float y_from_texture = from_start.texturePoint.y + i * y_step_size_from_texture;

        float x_to_texture = to_start.texturePoint.x + i * x_step_size_to_texture;
        float y_to_texture = to_start.texturePoint.y + i * y_step_size_to_texture;

        CanvasPoint from = CanvasPoint();
        CanvasPoint to = CanvasPoint();

        from.texturePoint = TexturePoint(x_from_texture, y_from_texture);
        to.texturePoint = TexturePoint(x_to_texture, y_to_texture);

        draw_texture_line(window,from, to, textureMap, affine_matrix);
    }
}

CanvasPoint find_mid_point(std::array<CanvasPoint, 3> vertices) {
    CanvasPoint mid_point;
    mid_point.y = vertices[1].y;

    auto top_bottom_x_diff = float(abs(int(vertices[0].x - vertices[2].x)));
    auto top_bottom_y_diff = float (abs(int(vertices[0].y - vertices[2].y)));

    float mid_x_diff = top_bottom_x_diff/top_bottom_y_diff * float(abs(int(vertices[0].y - vertices[1].y)));
    mid_point.depth = interpolateSingleFloats(vertices[0].depth, vertices[2].depth, 3)[1];
    if (vertices[0].x >= vertices[2].x) mid_point.x = vertices[0].x - mid_x_diff;
    else mid_point.x = vertices[0].x + mid_x_diff;
    return mid_point;
}

void fill_triangle(DrawingWindow &window, CanvasTriangle triangle, const Colour& colour) {
    std::array<CanvasPoint, 3> vertices = triangle.vertices;
    while (true) {
        if (vertices[0].y <= vertices[1].y && vertices[1].y <= vertices[2].y) break;
        else if (vertices[2].y < vertices[0].y) std::swap(vertices[2], vertices[0]);
        else if (vertices[1].y < vertices[0].y) std::swap(vertices[0], vertices[1]);
        else if (vertices[2].y < vertices[1].y) std::swap(vertices[2], vertices[1]);
    }

    CanvasPoint mid_point = find_mid_point(vertices);
    // draw top triangle
    fill_half_triangle(window, vertices[0], mid_point, vertices[1], colour);
    // draw bottom triangle
    fill_half_triangle(window, vertices[2], mid_point, vertices[1],colour);
    draw_line(window, vertices[1], mid_point, colour);
    draw_line(window, triangle.v0(), triangle.v1(), colour);
}

void draw_stroked_triangles(DrawingWindow &window, CanvasTriangle triangle, const Colour& colour) {
    draw_line(window, triangle.v0(), triangle.v1(), colour);
    draw_line(window, triangle.v1(), triangle.v2(), colour);
    draw_line(window, triangle.v2(), triangle.v0(), colour);
}

void draw_filled_triangles(DrawingWindow &window, CanvasTriangle triangle, const Colour& colour
                           ) {
    draw_stroked_triangles(window, triangle, colour);
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

CanvasPoint getCanvasIntersectionPoint(glm::vec3 cameraPosition,
                                       glm::vec3 vertexPosition,
                                       float focalLength,
                                       float scaling) {
    float camera_vertex_length = abs(cameraPosition.z - vertexPosition.z);
    float u = scaling*focalLength * (vertexPosition.x-cameraPosition.x) / camera_vertex_length;
    float v = -1*scaling*focalLength * (vertexPosition.y-cameraPosition.y) / camera_vertex_length;
    return {u+float(WIDTH)/2, v+float(HEIGHT)/2, vertexPosition.z};
}

void wire_frame_render(DrawingWindow &window,
                        const std::vector<ModelTriangle>& model_triangles,
                        glm::vec3 cameraPosition,
                        float focalLength,
                        float image_plane_scale) {

    for (auto & x : depth_buffer) {
        for (float & y : x) {
            y = 0.0;
        }
    }

    for (const auto& triangle : model_triangles) {

        std::vector<CanvasPoint> image_plane_triangle_vertices;
        for (auto vertex : triangle.vertices) {
            CanvasPoint sdl_vertex = getCanvasIntersectionPoint(cameraPosition, vertex, focalLength, image_plane_scale);
            image_plane_triangle_vertices.push_back(sdl_vertex);
        }
        CanvasTriangle image_plane_triangle = CanvasTriangle(image_plane_triangle_vertices[0],
                                                             image_plane_triangle_vertices[1],
                                                             image_plane_triangle_vertices[2]);

        draw_filled_triangles(window, image_plane_triangle,
                              triangle.colour);
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

glm::mat3x3 reverse_mtx(glm::mat3x3 mat) {
    glm::mat3x3 reverse_matrix;
    float determinant = 0;
    determinant += (mat[0][0] * (mat[1][(0 + 1) % 3] * mat[2][(0 + 2) % 3] - mat[1][(0 + 2) % 3] * mat[2][(0 + 1) % 3]));
    determinant += (mat[0][1] * (mat[1][(1 + 1) % 3] * mat[2][(1 + 2) % 3] - mat[1][(1 + 2) % 3] * mat[2][(1 + 1) % 3]));
    determinant += (mat[0][2] * (mat[1][(2 + 1) % 3] * mat[2][(2 + 2) % 3] - mat[1][(2 + 2) % 3] * mat[2][(2 + 1) % 3]));

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

    glm::mat3x3 result = reverse_mtx(texture_matrix) * sdl_matrix;

    return result;
}

void textureMapper(DrawingWindow &window, CanvasTriangle canvasTriangle, TextureMap textureMap) {
    std::array<CanvasPoint, 3> vertices = canvasTriangle.vertices;
    while (true) {
        if (vertices[0].y <= vertices[1].y && vertices[1].y <= vertices[2].y) break;
        else if (vertices[2].y <= vertices[0].y) std::swap(vertices[2], vertices[0]);
        else if (vertices[1].y <= vertices[0].y) std::swap(vertices[0], vertices[1]);
        else if (vertices[2].y <= vertices[1].y) std::swap(vertices[2], vertices[1]);
    }


    CanvasPoint midpoint = find_mid_point(vertices);
    std::array<CanvasPoint, 3> texture_vertices{
        CanvasPoint(vertices[0].texturePoint.x, vertices[0].texturePoint.y),
        CanvasPoint(vertices[1].texturePoint.x, vertices[1].texturePoint.y),
        CanvasPoint(vertices[2].texturePoint.x, vertices[2].texturePoint.y)
    };
    CanvasPoint mid_texture_point = find_mid_point(texture_vertices);


    midpoint.texturePoint = TexturePoint(mid_texture_point.x, mid_texture_point.y);

    glm::mat3x3 affineMtx = calculate_affine_mtx(canvasTriangle);
    fill_half_texture_triangle(window, vertices[0], vertices[0], midpoint, vertices[1], textureMap, affineMtx);
    fill_half_texture_triangle(window,  midpoint, vertices[1], vertices[2], vertices[2],textureMap, affineMtx);
    draw_stroked_triangles(window, canvasTriangle, Colour(255, 255, 255));
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
    std::vector<ModelTriangle> model_triangles = read_OBJ_files("cornell-box.obj", "cornell-box.mtl", 0.35);
    glm::vec3 camera_position = glm::vec3(0.0, 0.0, 4.0);
    float focal_length = 2.0;

    wire_frame_render(window, model_triangles, camera_position, focal_length, WIDTH / 2);
   /* Colour colour = Colour(rand()%255+1, rand()%255+1, rand()%255+1);
    CanvasTriangle triangle = CanvasTriangle(
            CanvasPoint(float(rand()%WIDTH+1), float(rand()%HEIGHT+1)),
            CanvasPoint(float(rand()%WIDTH+1), float(rand()%HEIGHT+1)),
            CanvasPoint(float(rand()%WIDTH+1), float(rand()%HEIGHT+1)));
    draw_filled_triangles(window, triangle, colour);*/
	while (true) {
		if (window.pollForInputEvents(event)) handleEvent(event, window);
        window.renderFrame();
	}
}
