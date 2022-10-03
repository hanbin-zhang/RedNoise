#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include "glm/vec3.hpp"

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
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
