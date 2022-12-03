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
#include "RayTriangleIntersection.h"
#include "ModelTriangle.h"
#include "glm/gtx/string_cast.hpp"

#define WIDTH 1280
#define HEIGHT 960

//glm::mat3 camera_orientation;
float depth_buffer[WIDTH][HEIGHT];
Colour colour_buffer[WIDTH][HEIGHT];
enum shadingType {Flat, Gourand, Phong, Nicht};
int shading = Flat;
bool isSoftShadow = false;
std::vector<glm::vec3> vertices;
glm::vec3 sphereCentre;
// corresponding texture file name
std::map<std::string, TextureMap> textureFilename;
std::vector<glm::vec3> thisLightCluster;

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

glm::mat3 calculate_affine_mtx(const ModelTriangle& triangle, int textureWidth, int  textureHeight) {
    glm::mat3 sdl_matrix = glm::mat3(
            glm::vec3 {triangle.vertices[0].x, triangle.vertices[0].y, triangle.vertices[0].z},
            glm::vec3 {triangle.vertices[1].x, triangle.vertices[1].y, triangle.vertices[1].z},
            glm::vec3 {triangle.vertices[2].x, triangle.vertices[2].y, triangle.vertices[2].z});
    glm::mat3 texture_matrix = glm::mat3(
            glm::vec3 {round(triangle.texturePoints[0].x*textureWidth), round(triangle.texturePoints[0].y*textureHeight), 1.0f},
            glm::vec3 {round(triangle.texturePoints[1].x*textureWidth), round(triangle.texturePoints[1].y*textureHeight), 1.0f},
            glm::vec3 {round(triangle.texturePoints[2].x*textureWidth), round(triangle.texturePoints[2].y*textureHeight), 1.0f});

    glm::mat3 result = texture_matrix * glm::inverse(sdl_matrix);
    return result;
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
        } else if (current_line.compare(0, 6, "map_Kd")==0) {
            textureFilename[current_colour_string] = TextureMap("../" + split(current_line, ' ')[1]);

        }
    }
    return colour_map;
}

std::vector<ModelTriangle> read_OBJ_files(const std::string& file_name,
                                          const std::string& colour_file_name,
                                          float scaling,
                                          glm::vec3 vertexShift) {

    std::string current_line;
    std::ifstream MyReadFile(file_name);
    std::vector<TexturePoint> textureVertices;
    vertices = std::vector<glm::vec3>();
    std::vector<ModelTriangle> triangles;
    Colour current_colour;
    std::map<std::string, Colour> colour_palette = read_colour_palette(colour_file_name);

    while (getline (MyReadFile, current_line)) {
        if (current_line.compare(0, 2, "v ")==0) {
            //put a vertex on to the list
            std::vector<std::string> vertices_string = split(current_line, ' ');
            vertices.emplace_back(std::stof(vertices_string[1]) * scaling + vertexShift[0],
                                  std::stof(vertices_string[2]) * scaling + vertexShift[1],
                                  std::stof(vertices_string[3]) * scaling + vertexShift[2]);

        } else if (current_line.compare(0, 2, "vt")==0) {
            std::vector<std::string> vertices_string = split(current_line, ' ');
            textureVertices.emplace_back(TexturePoint(std::stof(vertices_string[1]),
                                              std::stof(vertices_string[2])));
        } else if (current_line.compare(0, 6, "usemtl")==0) {
            current_colour = colour_palette[split(current_line, ' ')[1]];
        } else if (current_line.compare(0, 1, "f")==0) {
            // add a facet triangle
            std::vector<std::string> facets_string = split(current_line, ' ');
            ModelTriangle current_triangle;
            //`current_triangle.colour = Colour(255, 0, 0);
            current_triangle.colour = current_colour;
            for (int i = 1; i < int(facets_string.size()); ++i) {
                std::vector<std::string> facets_vertex_string = split(facets_string[i], '/');

                current_triangle.vertices[i-1] = vertices[std::stoi(facets_vertex_string[0])-1];

                if (textureFilename.count(current_colour.name)) {
                    current_triangle.texturePoints[i-1] = textureVertices[std::stoi(facets_vertex_string[1])-1];
                }
            }
            current_triangle.normal = glm::normalize(glm::cross((current_triangle.vertices[1] - current_triangle.vertices[0]),
                                                 (current_triangle.vertices[2] - current_triangle.vertices[0])));

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

float calculate_det(CanvasTriangle m_tria) {
    return (m_tria.v1().y - m_tria.v2().y) * (m_tria.v0().x - m_tria.v2().x) + (m_tria.v2().x - m_tria.v1().x) * (m_tria.v0().y - m_tria.v2().y);
}

glm::vec3 barycentricParams(CanvasTriangle m_tria, float x, float y) {
    float det_T = calculate_det(m_tria);
    float lam1 = 1/det_T * ((m_tria.v1().y - m_tria.v2().y) * (x-m_tria.v2().x) + (m_tria.v2().x - m_tria.v1().x) * (y - m_tria.v2().y));
    float lam2 = 1/det_T * ((m_tria.v2().y - m_tria.v0().y) * (x-m_tria.v2().x) + (m_tria.v0().x - m_tria.v2().x) * (y - m_tria.v2().y));
    float lam3 = 1 - lam1 - lam2;
    return glm::vec3 {lam1, lam2, lam3};
}

void draw_pixel(DrawingWindow &window, float x, float y,  const Colour& colour, CanvasTriangle m_tria, float det_T) {
    int round_x = int(round(x));
    int round_y = int(round(y));

    glm::vec3 lambdas = barycentricParams(m_tria, x, y);

    float lam1 = lambdas[0];
    float lam2 = lambdas[1];
    float lam3 = lambdas[2];

    float depth = lam1 * m_tria.v0().depth + lam2 * m_tria.v1().depth + lam3 * m_tria.v2().depth;

    if (!(round_x >= WIDTH || round_y >= HEIGHT || round_x <= 0 || round_y <= 0)) {
        if (depth>depth_buffer[round_x][round_y]) {
            if ( colour.green==0 && colour.red == 0 && colour.blue == 255) {
            }

            depth_buffer[round_x][round_y] = depth;
            colour_buffer[round_x][round_y] = colour;
            window.setPixelColour(size_t(round_x), size_t(round_y), colour_uint32(colour));
            }
        }
}

void draw_line_with_depth(DrawingWindow &window, CanvasPoint from, CanvasPoint to, const Colour& colour,
                          CanvasTriangle mother_triangle) {

    float det_T = calculate_det(mother_triangle);

    if (from.x == to.x && from.y == to.y) {
        draw_pixel(window, from.x, from.y, colour, mother_triangle, det_T);
        return;
    }

    float x_diff = to.x - from.x;

    float numberOfSteps = abs(x_diff);
    float x_step_size = x_diff / numberOfSteps;

    for (float i = 0; i <= (numberOfSteps); ++i) {

        float x = from.x + i*x_step_size;
        float y = from.y;

        draw_pixel(window, x, y,  colour, mother_triangle, det_T);

    }
}

void fill_half_triangle(DrawingWindow &window, CanvasPoint start,
                        CanvasPoint from_end, CanvasPoint to_end, const Colour& colour,
                        CanvasTriangle mother_triangle) {
    if ((from_end.x > WIDTH || from_end.x < 0) && (from_end.y > HEIGHT || from_end.y < 0) &&
    (to_end.x > WIDTH || to_end.x < 0) && (to_end.y > HEIGHT || to_end.y < 0) &&
    (start.x > WIDTH || start.x < 0) && (start.y > HEIGHT || start.y < 0))
        return;

    float x_diff_to = from_end.x - start.x;
    float x_diff_from = to_end.x - start.x;
    float y_diff = to_end.y - start.y;

    float numberOfSteps_x = std::max(abs(x_diff_from), abs(x_diff_to));
    float numberOfSteps = std::max(numberOfSteps_x, abs(y_diff));

    float x_step_size_from = x_diff_from / numberOfSteps;
    float x_step_size_to = x_diff_to / numberOfSteps;
    float y_step_size = y_diff / numberOfSteps;

    for (float i = 0.0; i <= numberOfSteps; ++i) {
        float x_from = start.x + i * x_step_size_from;
        float y_from = start.y + i * y_step_size;

        float x_to = start.x + i * x_step_size_to;
        float y_to = start.y + i * y_step_size;

        draw_line_with_depth(window, CanvasPoint(x_from, y_from),
                             CanvasPoint(x_to, y_to),
                             colour,
                             mother_triangle);

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

    auto top_bottom_x_diff = float(abs(vertices[0].x - vertices[2].x));
    auto top_bottom_y_diff = float (abs(vertices[0].y - vertices[2].y));
    auto top_bottom_depth_diff = float (abs(vertices[0].depth - vertices[2].depth));

    float mid_x_diff = top_bottom_x_diff/top_bottom_y_diff * float(abs(vertices[0].y - vertices[1].y));
    float mid_depth_diff = top_bottom_depth_diff/top_bottom_y_diff * float(abs(vertices[0].y - vertices[1].y));

    if (vertices[0].depth >= vertices[2].depth) mid_point.depth = vertices[0].depth-mid_depth_diff;
    else mid_point.depth = vertices[0].depth+mid_depth_diff;

    if (vertices[0].x >= vertices[2].x) mid_point.x = vertices[0].x - mid_x_diff;
    else mid_point.x = vertices[0].x + mid_x_diff;

    return mid_point;
}

void fill_triangle(DrawingWindow &window, CanvasTriangle triangle, const Colour& colour) {
    std::array<CanvasPoint, 3> vertices = triangle.vertices;
    if (vertices[2].y < vertices[0].y) std::swap(vertices[2], vertices[0]);
    if (vertices[1].y < vertices[0].y) std::swap(vertices[0], vertices[1]);
    if (vertices[2].y < vertices[1].y) std::swap(vertices[2], vertices[1]);
    CanvasPoint mid_point = find_mid_point(vertices);

    // draw top triangle

    fill_half_triangle(window, vertices[0], mid_point, vertices[1], colour, triangle);
    // draw bottom triangle
    fill_half_triangle(window, vertices[2], mid_point, vertices[1], colour, triangle);

}

void draw_stroked_triangles(DrawingWindow &window, CanvasTriangle triangle, const Colour& colour) {
    draw_line(window, triangle.v0(), triangle.v1(), colour);
    draw_line(window, triangle.v1(), triangle.v2(), colour);
    draw_line(window, triangle.v2(), triangle.v0(), colour);
}

void draw_filled_triangles(DrawingWindow &window, CanvasTriangle triangle, const Colour& colour
                           ) {
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

glm::mat3 lookAt(glm::vec3 cameraPosition) {
    glm::vec3 forward = glm::normalize( cameraPosition);
    glm::vec3 right = glm::normalize(glm::cross(glm::vec3(0, 1, 0), forward));
    glm::vec3 up = glm::normalize(glm::cross(forward, right));
    glm::mat3 camera_orbit_orientation = {right, up, forward};
    return camera_orbit_orientation;
}

bool closestIntersectionTests(glm::vec3 possibleSolution) {
    return (((possibleSolution[1] >= 0.0) && (possibleSolution[1] <= 1.0)) &&
            ((possibleSolution[2] >= 0.0) && (possibleSolution[2] <= 1.0)) &&
            ((possibleSolution[1] + possibleSolution[2]) <= 1.0) &&
            (possibleSolution[0] > 0.0001));

}

glm::vec3 calculateVertexNormal(const std::vector<ModelTriangle>& model_triangles,
                                glm::vec3 vertex) {
    int involveFaceNum = 0;
    glm::vec3 involveFaceSum = glm::vec3 {0, 0, 0};

    for (const auto& triangle : model_triangles) {
        if (triangle.vertices[0] == vertex ||
            triangle.vertices[1] == vertex ||
            triangle.vertices[2] == vertex ){
            involveFaceNum++;
            involveFaceSum += triangle.normal;
        }
    }
    return glm::normalize(involveFaceSum / float (involveFaceNum));
}

float proximityParameter(glm::vec3 lightSource, glm::vec3 vertexPosition, float lightIntensity) {
    float distance = glm::length(lightSource - vertexPosition);

    // 1/4πr2
    float para = lightIntensity / float (4 * M_PI * distance * distance);
    if (para <= 0.0) return 0.0;
    else if (para > 1.0) return 1.0;
    else return para;
}

float angleOfIncidentParam(glm::vec3 lightSource, glm::vec3 vertex,glm::vec3 targetNormal) {
    glm::vec3 toLightDirection = glm::normalize(lightSource - vertex);
    return glm::clamp<float>(
            glm::dot(targetNormal, toLightDirection),
            0.0, 1.0);
}

glm::vec3 calculateReflection(glm::vec3 ri, glm::vec3 targetNormal) {
    return ri - float (2.0) * targetNormal * (glm::dot(ri, targetNormal));
}

float specularParam(glm::vec3 lightSource, glm::vec3 cameraPosition, glm::vec3 targetNormal, glm::vec3 vertex) {
    glm::vec3 ri = glm::normalize(vertex - lightSource);
    glm::vec3 reflection = calculateReflection(ri, targetNormal);
    glm::vec3 view = glm::normalize(cameraPosition - vertex);
    auto param = glm::clamp<float>(glm::dot(reflection, view), 0.0, 1.0);
    return param;

}

float lightParam(glm::vec3 lightSource, glm::vec3 cameraPosition, glm::vec3 vertex,
                 const std::vector<ModelTriangle>& model_triangles, glm::vec3 targetNormal) {

    float proximityParam = proximityParameter(lightSource,
                                              vertex, 20.0);
    float aoIParam = angleOfIncidentParam( lightSource, vertex, targetNormal);
    float specularP = specularParam( lightSource, cameraPosition, targetNormal, vertex);
    return glm::clamp<float>( proximityParam * aoIParam + float (pow(specularP, 512)) + 0.2f, 0.0, 1.0);
}

float gouraudLight(const std::vector<ModelTriangle>& model_triangles,
                   const RayTriangleIntersection& intersection,
                   glm::vec3 lightSource, glm::vec3 cameraPosition) {

    glm::vec3 vertex = intersection.intersectionPoint;
    glm::vec3 lambdas;
    if (intersection.intersectedTriangle.vertices[0].x == intersection.intersectedTriangle.vertices[1].x &&
        intersection.intersectedTriangle.vertices[0].x == intersection.intersectedTriangle.vertices[2].x &&
            intersection.intersectedTriangle.vertices[1].x == intersection.intersectedTriangle.vertices[2].x) {
        CanvasTriangle canvasTriangle = CanvasTriangle(CanvasPoint(intersection.intersectedTriangle.vertices[0].y,
                                                                   intersection.intersectedTriangle.vertices[0].z),
                                                       CanvasPoint(intersection.intersectedTriangle.vertices[1].y,
                                                                   intersection.intersectedTriangle.vertices[1].z),
                                                       CanvasPoint(intersection.intersectedTriangle.vertices[2].y,
                                                                   intersection.intersectedTriangle.vertices[2].z));
        lambdas = barycentricParams(canvasTriangle,vertex[1], vertex[2]);
    } else {
        CanvasTriangle canvasTriangle = CanvasTriangle(CanvasPoint(intersection.intersectedTriangle.vertices[0].x,
                                                                   intersection.intersectedTriangle.vertices[0].y),
                                                       CanvasPoint(intersection.intersectedTriangle.vertices[1].x,
                                                                   intersection.intersectedTriangle.vertices[1].y),
                                                       CanvasPoint(intersection.intersectedTriangle.vertices[2].x,
                                                                   intersection.intersectedTriangle.vertices[2].y));
        lambdas = barycentricParams(canvasTriangle,vertex[0], vertex[1]);

    }

    glm::vec3 normalV0 = calculateVertexNormal(model_triangles, intersection.intersectedTriangle.vertices[0]);
    glm::vec3 normalV1 = calculateVertexNormal(model_triangles, intersection.intersectedTriangle.vertices[1]);
    glm::vec3 normalV2 = calculateVertexNormal(model_triangles, intersection.intersectedTriangle.vertices[2]);


    float lightV0 = lightParam(lightSource, cameraPosition, intersection.intersectedTriangle.vertices[0],
                               model_triangles, normalV0);
    float lightV1 = lightParam(lightSource, cameraPosition, intersection.intersectedTriangle.vertices[1],
                               model_triangles, normalV1);
    float lightV2 = lightParam(lightSource, cameraPosition, intersection.intersectedTriangle.vertices[2],
                               model_triangles, normalV2);

    return glm::clamp<float>((lightV0*lambdas[0] + lightV1*lambdas[1] + lightV2*lambdas[2]), 0.0, 1.0);
}

glm::vec3 eachVertexNormal(RayTriangleIntersection intersection, glm::vec3 vertex,
                           const std::vector<ModelTriangle>& model_triangles) {
    glm::vec3 lambdas;
    if (intersection.intersectedTriangle.vertices[0].x == intersection.intersectedTriangle.vertices[1].x &&
        intersection.intersectedTriangle.vertices[0].x == intersection.intersectedTriangle.vertices[2].x &&
        intersection.intersectedTriangle.vertices[1].x == intersection.intersectedTriangle.vertices[2].x) {
        CanvasTriangle canvasTriangle = CanvasTriangle(CanvasPoint(intersection.intersectedTriangle.vertices[0].y,
                                                                   intersection.intersectedTriangle.vertices[0].z),
                                                       CanvasPoint(intersection.intersectedTriangle.vertices[1].y,
                                                                   intersection.intersectedTriangle.vertices[1].z),
                                                       CanvasPoint(intersection.intersectedTriangle.vertices[2].y,
                                                                   intersection.intersectedTriangle.vertices[2].z));
        lambdas = barycentricParams(canvasTriangle,vertex[1], vertex[2]);
    } else {
        CanvasTriangle canvasTriangle = CanvasTriangle(CanvasPoint(intersection.intersectedTriangle.vertices[0].x,
                                                                   intersection.intersectedTriangle.vertices[0].y),
                                                       CanvasPoint(intersection.intersectedTriangle.vertices[1].x,
                                                                   intersection.intersectedTriangle.vertices[1].y),
                                                       CanvasPoint(intersection.intersectedTriangle.vertices[2].x,
                                                                   intersection.intersectedTriangle.vertices[2].y));
        lambdas = barycentricParams(canvasTriangle,vertex[0], vertex[1]);

    }

    glm::vec3 normalV0 = calculateVertexNormal(model_triangles, intersection.intersectedTriangle.vertices[0]);
    glm::vec3 normalV1 = calculateVertexNormal(model_triangles, intersection.intersectedTriangle.vertices[1]);
    glm::vec3 normalV2 = calculateVertexNormal(model_triangles, intersection.intersectedTriangle.vertices[2]);

    return lambdas[0]*normalV0 + lambdas[1]*normalV1 + lambdas[2]*normalV2;
}

float phongLight(const std::vector<ModelTriangle>& model_triangles,
                   const RayTriangleIntersection& intersection,
                   glm::vec3 lightSource, glm::vec3 cameraPosition) {

    glm::vec3 vertex = intersection.intersectionPoint;

    glm::vec3 targetNormal = eachVertexNormal(intersection, vertex, model_triangles);

    return lightParam(lightSource, cameraPosition, vertex, model_triangles, targetNormal);
}

glm::vec3 mirror(glm::vec3 vertex, const std::vector<ModelTriangle>& model_triangles, glm::vec3 cameraPosition,
              const RayTriangleIntersection& intersection) {
    glm::vec3 reflection = glm::normalize(calculateReflection(glm::normalize(vertex - cameraPosition), intersection.intersectedTriangle.normal));

    return reflection;
}

glm::vec3 refract(glm::vec3 incident, glm::vec3 normal, float refractive_index) {
    float incidentAngle = glm::dot(incident, normal);
    float n1 = 1, n2 = refractive_index;
    float ratio = n1 / n2;
    if (incidentAngle > 0) {
        ratio = 1/ratio;
        normal = -normal;
    }
    incidentAngle = abs(incidentAngle);
    float  k = 1 - ratio * ratio * (1 - incidentAngle * incidentAngle);
    return k < 0 ? glm::vec3 {0, 0, 0} :glm::normalize(incident * ratio + normal * (ratio * incidentAngle - sqrtf(k)));
}

float fresnelLaw(glm::vec3 incident, glm::vec3 normal, float refractiveIndex) {

    float cosi = glm::dot(incident, normal);
    float etai = 1;
    float etat = refractiveIndex;

    if (cosi > 0) {

        std::swap(etai, etat);
    }
    float eta = etai / etat;
    float sint = eta * sqrtf(std::max(0.f, 1 - cosi * cosi));
    if (sint >= 1) return 1;
    else {
        float cost = sqrtf(std::max(0.f, 1 - sint * sint));
        cosi = fabsf(cosi);

        float rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        float rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        return (rs * rs + rp * rp )*0.5f;
    }
}

RayTriangleIntersection getClosestIntersection(glm::vec3 camera_position, glm::vec3 ray_direction,
                                               const std::vector<ModelTriangle>& triangles) {
    RayTriangleIntersection intersection;

    auto absolute_distance = float (INT32_MAX-1);
    for (int i = 0; i < int (triangles.size()); ++i) {
        glm::vec3 e0 = triangles[i].vertices[1] - triangles[i].vertices[0];
        glm::vec3 e1 = triangles[i].vertices[2] - triangles[i].vertices[0];
        glm::vec3 SPVector = camera_position - triangles[i].vertices[0];
        glm::mat3 DEMatrix(-ray_direction, e0, e1);
        glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

        if (closestIntersectionTests(possibleSolution)) {
            if (possibleSolution[0] < absolute_distance) {
                absolute_distance = possibleSolution[0];
                glm::vec3 r = camera_position + ray_direction * possibleSolution[0];
                intersection = RayTriangleIntersection(r, possibleSolution[0], triangles[i], std::size_t(i));
            }
        }
    }
    return intersection;
}

std::vector<glm::vec3> lightCluster(glm::vec3 lightSource, int radian, float stepSize) {
    std::vector<glm::vec3> lights;
    for (int i = -radian; i <= radian ; i++) {
        for (int j = -radian; j <= radian ; j++) {
            if ((i * i + j * j)<=radian*radian) {
                glm::vec3 oneLight = {lightSource.x+float (i)*stepSize,
                                      lightSource.y,
                                      lightSource.z+float (j)*stepSize};

                lights.emplace_back(oneLight);
            }
        }
    }
    return lights;
}

float softShadowParam(glm::vec3 point,
                      const std::vector<ModelTriangle>& triangles,
                      const RayTriangleIntersection& intersection) {
    int lightNumber = 0;
    for (int i = 0; i < int (thisLightCluster.size()); ++i) {
        glm::vec3 light = thisLightCluster[i];
        glm::vec3 fromLightDirection = glm::normalize(point - light);

        RayTriangleIntersection lightIntersection =
                getClosestIntersection(light, fromLightDirection, triangles);

        if (intersection.triangleIndex == lightIntersection.triangleIndex) {
            lightNumber = lightNumber + 1;
        } else if (lightIntersection.intersectedTriangle.colour.name.compare(0, 3, "Red") == 0) {
            return 0.7;
        }
    }
        return float (lightNumber) / float (thisLightCluster.size());
}

Colour shootRay(glm::vec3 cameraPosition,
                glm::vec3 rayDirection,
                glm::vec3 lightSource,
                std::vector<ModelTriangle> triangles,
                int recurrentNumber) {
    if (recurrentNumber >= 5) {
        return {255, 255, 255};
    }
    RayTriangleIntersection intersection = getClosestIntersection(cameraPosition, rayDirection, triangles);
    Colour targetColour = intersection.intersectedTriangle.colour;
    if (textureFilename.count(intersection.intersectedTriangle.colour.name)) {
        TextureMap textureMap = textureFilename[intersection.intersectedTriangle.colour.name];
        glm::mat3 affineMat = calculate_affine_mtx(intersection.intersectedTriangle, float (textureMap.width), float (textureMap.height));
        glm::vec3 texturePoint = affineMat * intersection.intersectionPoint;

        uint32_t colourInt = textureMap.pixels[round(texturePoint[0]) +
                                               round(texturePoint[1]) * textureMap.width];
        unsigned  mask;
        mask = 0xff;
        uint32_t red = (colourInt >> 16) & mask;
        uint32_t green = (colourInt >> 8) & mask;
        uint32_t blue = (colourInt) & mask;
        targetColour = Colour {int(red), int(green), int (blue)};
    }

    if (intersection.intersectedTriangle.colour.name.compare(0, 5, "Green")==0) {
        TextureMap textureMap = textureFilename["mars"];
        glm::vec3 spherePoint = (intersection.intersectionPoint - sphereCentre) / 0.35f;
        float faiz = atan2(spherePoint.z, spherePoint.x);
        float theta = asin(spherePoint.y);

        float u = 1.0f - (faiz + M_PI) / (2.0f*M_PI);
        float v = (theta + M_PI/2.0f) / M_PI;

        int x = int (u * textureMap.width);
        int y = int ((1-v) * textureMap.height);

        uint32_t colourInt = textureMap.pixels[x +
                                               y * textureMap.width];
        unsigned  mask;
        mask = 0xff;
        uint32_t red = (colourInt >> 16) & mask;
        uint32_t green = (colourInt >> 8) & mask;
        uint32_t blue = (colourInt) & mask;
        targetColour = Colour {int(red), int(green), int (blue)};
    }

    if (intersection.intersectedTriangle.colour.name.compare(0, 7, "Magenta")==0) {
        glm::vec3 reflection = mirror(intersection.intersectionPoint,
                                      triangles,
                                      cameraPosition,
                                      intersection);

        return shootRay(intersection.intersectionPoint,
                                reflection,
                                lightSource,
                                triangles,
                                recurrentNumber+1);
    } else if (intersection.intersectedTriangle.colour.name.compare(0, 3, "Red")==0) {


        glm::vec3 reflection = mirror(intersection.intersectionPoint,
                                      triangles,
                                      cameraPosition,
                                      intersection);

        Colour reflectionColour = shootRay(intersection.intersectionPoint,
                                           reflection,
                                           lightSource,
                                           triangles,
                                           recurrentNumber+1);

        intersection.intersectedTriangle.normal = intersection.intersectedTriangle.normal*-1.0f;

        float refractiveIndex = 1.3;
        glm::vec3 refractDirection = refract(rayDirection,
                                             intersection.intersectedTriangle.normal,
                                             refractiveIndex);

        if (refractDirection == glm::vec3(0, 0, 0)) {
            return reflectionColour;
        }

        float incidentAngle = glm::dot(rayDirection, intersection.intersectedTriangle.normal);
        glm::vec3 updatePoint;

        if (incidentAngle < 0) {
            updatePoint = intersection.intersectionPoint - intersection.intersectedTriangle.normal * 0.0001f;
        } else {
            updatePoint = intersection.intersectionPoint + intersection.intersectedTriangle.normal * 0.0001f;
        }

        Colour refraColour = shootRay(updatePoint,
                                      refractDirection,
                                      lightSource,
                                      triangles,
                                      recurrentNumber+1);

        float reflectiveConst = fresnelLaw(rayDirection,
                                           intersection.intersectedTriangle.normal, refractiveIndex);
        float refractiveConst = 1 - reflectiveConst;

        Colour doubleRcolour;
        doubleRcolour.red = (reflectiveConst * reflectionColour.red) + (refractiveConst * refraColour.red);
        doubleRcolour.green = (reflectiveConst * reflectionColour.green) + (refractiveConst * refraColour.green);
        doubleRcolour.blue = (reflectiveConst * reflectionColour.blue) + (refractiveConst * refraColour.blue);
        return doubleRcolour;
    }
    // shadow

    float light_param;
    if (shading == Flat)
        light_param = lightParam(lightSource, cameraPosition, intersection.intersectionPoint,
                                 triangles, intersection.intersectedTriangle.normal);
    else if (shading == Gourand)
        light_param = gouraudLight(triangles, intersection,
                                   lightSource, cameraPosition);
    else if (shading == Phong)
        light_param = phongLight(triangles, intersection,
                                 lightSource, cameraPosition);
    else if (shading == Nicht) light_param = 1;

    if (isSoftShadow) {
        float ssParam = softShadowParam(intersection.intersectionPoint, triangles, intersection);
        ssParam = glm::clamp<float>(ssParam + 0.2, 0.0f, 1.0f);
        targetColour = Colour(targetColour.red * (light_param) * ssParam,
                              (targetColour.green) * (light_param) * ssParam,
                              (targetColour.blue) * (light_param)* ssParam);
    } else {
        glm::vec3 fromLightDirection = glm::normalize(intersection.intersectionPoint - lightSource);

        RayTriangleIntersection lightIntersection =
                getClosestIntersection(lightSource, fromLightDirection, triangles);

        if (intersection.triangleIndex == lightIntersection.triangleIndex) {

            targetColour = Colour(targetColour.red * (light_param),
                                  (targetColour.green) * (light_param),
                                  (targetColour.blue) * (light_param)
            );
        } else {
            float shadowPram;
            if (lightIntersection.intersectedTriangle.colour.name.compare(0, 3, "Red") == 0) shadowPram = 0.45;
            else shadowPram = 0.2;

            targetColour = Colour(float(targetColour.red) * shadowPram,
                                  float(targetColour.green) * shadowPram,
                                  float(targetColour.blue) * shadowPram
            );
        }
    }
    return targetColour;
}


void rayTracingRender(DrawingWindow &window,
                      const std::vector<ModelTriangle>& model_triangles,
                      glm::vec3 cameraPosition,
                      float focalLength,
                      float scaling,
                      float orbiting_radian,
                      glm::vec3 lightSource
                      ) {
    glm::mat3 orbiting =  {cos(orbiting_radian), 0, sin(orbiting_radian),
                           0, 1, 0,
                           -sin(orbiting_radian), 0, cos(orbiting_radian)};

    cameraPosition = cameraPosition * orbiting;

    glm::mat3 camera_orbit_orientation = lookAt(cameraPosition);

    for (int u = 0; u < WIDTH; ++u) {
        float x = (float(u) - float(WIDTH) / 2) / scaling ;
        for (int v = 0; v < HEIGHT; ++v) {
            float y = -1 * (float(v) - float(HEIGHT) / 2) / scaling ;
            glm::vec3 image_plane_vertex = glm::vec3{x + cameraPosition.x, y + cameraPosition.y,
                                                     cameraPosition.z - focalLength};
            glm::vec3 imagePlaneDirection = glm::normalize(image_plane_vertex - cameraPosition);

            imagePlaneDirection = glm::normalize(imagePlaneDirection  * glm::inverse(camera_orbit_orientation) );

            Colour colour = shootRay(cameraPosition,
                                     imagePlaneDirection,
                                     lightSource,
                                     model_triangles,
                                     0);

            window.setPixelColour(std::size_t(u), std::size_t(v),
                                  colour_uint32(colour));

        }
    }
}

CanvasPoint getCanvasIntersectionPoint(glm::vec3 cameraPosition,
                                       glm::vec3 vertexPosition,
                                       float focalLength,
                                       float scaling,
                                       glm::mat3 camera_orbit_orientation) {

    glm::vec3 diff =  vertexPosition - cameraPosition;
    diff = diff * camera_orbit_orientation;

    float u = float (WIDTH) - round(float(WIDTH)/2 + scaling*focalLength * (diff.x) / diff.z);
    float v = round(float(HEIGHT)/2 + scaling*focalLength * (diff.y) / diff.z);
    return {u, v, abs(1/diff.z)} ;
}

void initialize_depth_buffer() {
    for (auto & x : depth_buffer) {
        for (float & y : x) {
            y = 0;
        }
    }
}

void Rasterised_render(DrawingWindow &window,
                       const std::vector<ModelTriangle>& model_triangles,
                       glm::vec3 cameraPosition,
                       float focalLength,
                       float image_plane_scale,
                       float orbiting_radian,
                       bool is_wire_frame) {

    initialize_depth_buffer();

    glm::mat3 orbiting =  {cos(orbiting_radian), 0, sin(orbiting_radian),
                           0, 1, 0,
                           -sin(orbiting_radian), 0, cos(orbiting_radian)};

    cameraPosition = orbiting*cameraPosition;

    glm::mat3 camera_orbit_orientation = lookAt(cameraPosition);

    for (const auto& triangle : model_triangles) {

        std::vector<CanvasPoint> image_plane_triangle_vertices;
        for (auto vertex : triangle.vertices) {

            CanvasPoint sdl_vertex = getCanvasIntersectionPoint(cameraPosition, vertex, focalLength, image_plane_scale,
                                                                camera_orbit_orientation);
            image_plane_triangle_vertices.push_back(sdl_vertex);
        }

        CanvasTriangle image_plane_triangle = CanvasTriangle(image_plane_triangle_vertices[0],
                                                             image_plane_triangle_vertices[1],
                                                             image_plane_triangle_vertices[2]);

        if (is_wire_frame) {
            draw_stroked_triangles(window, image_plane_triangle,
                                   triangle.colour);
        } else {
            draw_filled_triangles(window, image_plane_triangle,
                                  triangle.colour);
        }
    }
}

//void textureMapper(DrawingWindow &window, CanvasTriangle canvasTriangle, TextureMap textureMap) {
//    std::array<CanvasPoint, 3> vertices = canvasTriangle.vertices;
//
//     if (vertices[2].y < vertices[0].y) std::swap(vertices[2], vertices[0]);
//     if (vertices[1].y < vertices[0].y) std::swap(vertices[0], vertices[1]);
//     if (vertices[2].y < vertices[1].y) std::swap(vertices[2], vertices[1]);
//
//
//    CanvasPoint midpoint = find_mid_point(vertices);
//    std::array<CanvasPoint, 3> texture_vertices{
//        CanvasPoint(vertices[0].texturePoint.x, vertices[0].texturePoint.y),
//        CanvasPoint(vertices[1].texturePoint.x, vertices[1].texturePoint.y),
//        CanvasPoint(vertices[2].texturePoint.x, vertices[2].texturePoint.y)
//    };
//    CanvasPoint mid_texture_point = find_mid_point(texture_vertices);
//
//
//    midpoint.texturePoint = TexturePoint(mid_texture_point.x, mid_texture_point.y);
//
//    glm::mat3x3 affineMtx = calculate_affine_mtx(canvasTriangle);
//    fill_half_texture_triangle(window, vertices[0], vertices[0], midpoint, vertices[1], textureMap, affineMtx);
//    fill_half_texture_triangle(window,  midpoint, vertices[1], vertices[2], vertices[2],textureMap, affineMtx);
//    draw_stroked_triangles(window, canvasTriangle, Colour(255, 255, 255));
//}

void handleEvent(SDL_Event event, DrawingWindow &window, glm::vec3* camera_position,
                 float* x_rotate, float* y_rotate, bool* is_rotate, int* render_mode,
                 float* lightX, float* lightY, float* lightZ) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) camera_position->x += -0.1;
		else if (event.key.keysym.sym == SDLK_RIGHT) camera_position->x += 0.1;
		else if (event.key.keysym.sym == SDLK_UP) camera_position->y += -0.1;
		else if (event.key.keysym.sym == SDLK_DOWN) camera_position->y += 0.1;
        else if (event.key.keysym.sym == SDLK_w) *x_rotate += float (M_PI/180);
        else if (event.key.keysym.sym == SDLK_s) *x_rotate -= float (M_PI/180);
        else if (event.key.keysym.sym == SDLK_a) *y_rotate -= float (M_PI/180);
        else if (event.key.keysym.sym == SDLK_d) *y_rotate += float (M_PI/180);
        else if (event.key.keysym.sym == SDLK_p) *is_rotate = ! *is_rotate;
        else if (event.key.keysym.sym == SDLK_r) *render_mode = 0;
        else if (event.key.keysym.sym == SDLK_t) *render_mode = 1;
        else if (event.key.keysym.sym == SDLK_y) *render_mode = 2;
        else if (event.key.keysym.sym == SDLK_f) *lightX += 0.1;
        else if (event.key.keysym.sym == SDLK_v) *lightX -= 0.1;
        else if (event.key.keysym.sym == SDLK_g) *lightY += 0.1;
        else if (event.key.keysym.sym == SDLK_b) *lightY -= 0.1;
        else if (event.key.keysym.sym == SDLK_h) *lightZ += 0.1;
        else if (event.key.keysym.sym == SDLK_n) *lightZ -= 0.1;
        else if (event.key.keysym.sym == SDLK_c) isSoftShadow = !isSoftShadow;
        else if (event.key.keysym.sym == SDLK_1) shading = Flat;
        else if (event.key.keysym.sym == SDLK_2) shading = Gourand;
        else if (event.key.keysym.sym == SDLK_3) shading = Phong;
        else if (event.key.keysym.sym == SDLK_4) shading = Nicht;

    } else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

//void calculate_camera_orientation(float x_rotate_radian, float y_rotate_radian) {
//    camera_orientation = {glm::vec3 {1, 0, 0},
//                          glm::vec3 {0, 1, 0},
//                          glm::vec3 {0, 0, 1},};
//
//    glm::mat3 x_rotate_mat= {glm::vec3{1, 0, 0},
//                             glm::vec3{0, cos(x_rotate_radian), -sin(x_rotate_radian)},
//                             glm::vec3{0, sin(x_rotate_radian), cos(x_rotate_radian)}};
//
//    glm::mat3 y_rotate_mat= {glm::vec3{cos(y_rotate_radian), 0, sin(y_rotate_radian)},
//                             glm::vec3{0, 1, 0},
//                             glm::vec3{-sin(y_rotate_radian), 0, cos(x_rotate_radian)}};
//
//    glm::vec3 camera_position;
//
//    camera_orientation = x_rotate_mat * y_rotate_mat * camera_orientation;
//}

void calculateSphereCentre() {
    glm::vec3 v1 = vertices[0], v2 = vertices[10], v3 = vertices[20], v4 = vertices[30];

    glm::mat3 A = 2.0f * glm::mat3 {glm::vec3 {v2.x -v1.x, v3.x - v2.x, v4.x - v3.x},
                   glm::vec3 {v2.y - v1.y, v3.y - v2.y, v4.y - v3.y},
                   glm::vec3 {v2.z - v1.z, v3.z - v2.z, v4.z - v3.z}};

    glm::vec3 beta = -1.0f * glm::vec3 {v1.x*v1.x - v2.x*v2.x + v1.y*v1.y - v2.y*v2.y + v1.z*v1.z - v2.z*v2.z,
                                        v2.x*v2.x - v3.x*v3.x + v2.y*v2.y - v3.y*v3.y + v2.z*v2.z - v3.z*v3.z,
                                        v3.x*v3.x - v4.x*v4.x + v3.y*v3.y - v4.y*v4.y + v3.z*v3.z - v4.z*v4.z};
    sphereCentre = glm::inverse(A) * beta;
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
    std::vector<ModelTriangle> box_model_triangles = read_OBJ_files("../textured-cornell-box.obj",
                                                                "../textured-cornell-box.mtl",
                                                                0.35,
                                                                glm::vec3(0, 0, 0));
    std::vector<ModelTriangle> sphere_model_triangles = read_OBJ_files("../sphere.obj",
                                                                       "../textured-cornell-box.mtl",
                                                                       0.35,
                                                                       glm::vec3(-0.5, -1.1349, 0.0));
    textureFilename["mars"] = TextureMap("../mars.ppm");
    std::vector<ModelTriangle> model_triangles;
    for (const auto& triangle : box_model_triangles) {
        model_triangles.emplace_back(triangle);
    }
    for (const auto& triangle : sphere_model_triangles) {
        model_triangles.emplace_back(triangle);
    }
    calculateSphereCentre();
    std::cout << glm::to_string(sphereCentre) << std::endl;
    glm::vec3 initial_camera_position = glm::vec3(0.0, 0.0, 4.0);
    float focal_length = 2.0;
    float x_rotate_radian = 0;
    float y_rotate_radian = 0;
    float orbiting_radian = 0;
    bool is_rotate = false;
    int render_mode = 0;
    float lightX = 0.0;
    float lightY = 0.5;
    float lightZ = 0.3;
    glm::vec3 lightSource = {lightX, lightY, lightZ};
    thisLightCluster = lightCluster(lightSource, 5, 0.1);

    while (true) {
        lightSource = {lightX, lightY, lightZ};
		if (window.pollForInputEvents(event)) handleEvent(event, window,
                                                          &initial_camera_position, &x_rotate_radian, &y_rotate_radian,
                                                          &is_rotate, &render_mode ,&lightX, &lightY, &lightZ);
        window.clearPixels();

        //calculate_camera_orientation(x_rotate_radian, y_rotate_radian);
        if (is_rotate) {
            if (orbiting_radian >= M_PI*2) orbiting_radian = 0;
            orbiting_radian += M_PI / 144;
        }

        switch (render_mode) {
            case 1:
                Rasterised_render(window, model_triangles,
                                  initial_camera_position,
                                  focal_length, WIDTH / 2,
                                  orbiting_radian,
                                  true);
                break;
            case 2:
                rayTracingRender(window, model_triangles,
                                 initial_camera_position,
                                 focal_length, WIDTH / 2,
                                 orbiting_radian,
                                 lightSource);
                break;
            default:
                Rasterised_render(window, model_triangles,
                                  initial_camera_position,
                                  focal_length, WIDTH / 2,
                                  orbiting_radian,
                                  false);
                break;
        }

        window.renderFrame();
	}
}
