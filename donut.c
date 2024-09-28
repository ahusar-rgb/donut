#include <stdio.h>
#include <stdlib.h>
#include<unistd.h>
#include <math.h>
#include <time.h>

#define SIZE 75
// #define SIZE 150
#define CIRCLES 80
#define POINTS_IN_CIRCLE 80
#define PI 3.14159
#define XY_FACTOR 2

typedef struct tensor {
    char*** array;
    int x_size;
    int y_size;
    int z_size;
} Tensor;


char pixels[12] = {'@', '$', '#', '*', '!', '=', ';', ':', '~', '-', ',', '.'};
int lightVect[] = {1, 1, 1};


Tensor allocateTensor(int x_size, int y_size, int z_size) {
    char*** points = (char***)malloc(x_size * sizeof(char**));
    for (int i = 0; i < x_size; i++) {
        points[i] = (char**)malloc(y_size * sizeof(char*));
        for (int j = 0; j < y_size; j++) {
            points[i][j] = (char*)malloc(z_size * sizeof(char));
        }
    }
    Tensor result = {points, x_size, y_size, z_size};
    return result;
}

void fillTensor(Tensor tensor, char fill) {
    for (int i = 0; i < tensor.x_size; i++) {
        for (int j = 0; j < tensor.y_size; j++) {
            for (int k = 0; k < tensor.z_size; k++) {
                tensor.array[i][j][k] = fill;
            }
        }
    }
}

void freeTensor(Tensor tensor) {
    for (int i = 0; i < tensor.x_size; i++) {
        for (int j = 0; j < tensor.y_size; j++) {
            free(tensor.array[i][j]);
        }
        free(tensor.array[i]);
    }
    free(tensor.array);
}

char findPixel(Tensor tensor, int x, int y, int z) {
    if (tensor.array[x][y][z] != ' ') {
        return tensor.array[x][y][z];
    }

    int empty_count = 0;
    int counted = 0;
    if (x > 0) {
        empty_count += tensor.array[x - 1][y][z] == ' ' ? 1 : 0;
        counted++;
    }
    if (x < tensor.x_size - 1) {
        empty_count += tensor.array[x + 1][y][z] == ' ' ? 1 : 0;
        counted++;
    }
    if (y > 0) {
        empty_count += tensor.array[x][y - 1][z] == ' ' ? 1 : 0;
        counted++;
    }
    if (y < tensor.y_size - 1) {
        empty_count += tensor.array[x][y + 1][z] == ' ' ? 1 : 0;
        counted++;
    }
    if (x > 0 && y > 0) {
        empty_count += tensor.array[x - 1][y - 1][z] == ' ' ? 1 : 0;
        counted++;
    }
    if (x > 0 && y < tensor.y_size - 1) {
        empty_count += tensor.array[x - 1][y + 1][z] == ' ' ? 1 : 0;
        counted++;
    }
    if (x < tensor.x_size - 1 && y > 0) {
        empty_count += tensor.array[x + 1][y - 1][z] == ' ' ? 1 : 0;
        counted++;
    }
    if (x < tensor.x_size - 1 && y < tensor.y_size - 1) {
        empty_count += tensor.array[x + 1][y + 1][z] == ' ' ? 1 : 0;
        counted++;
    }

    if (counted - empty_count >= 2) {
        return '.';
    }
    else {
        return ' ';
    }
}

void printTensorProj(Tensor tensor) {
    for (int x = 0; x < tensor.x_size; x += 2) {
        for(int y = 0; y < tensor.y_size; y++) {
            for(int z = 0; z < tensor.z_size; z++) {
                // char pixel = findPixel(tensor, x, y, z);
                char pixel = tensor.array[x][y][z];
                if (pixel != ' ') {
                    printf("%c", pixel);
                    break;
                }
        
                if (z == tensor.z_size - 1) {
                    printf(" ");
                }   
            }
        }
        printf("\n");
    }
}

void clearScreen() {
    printf("\033[H\033[J");
}

int* rotatePoint(int* point, int* pivot, char axis, double a) {
    int x = point[0] - pivot[0];
    int y = point[1] - pivot[1];
    int z = point[2] - pivot[2];

    int* result = (int*)malloc(3 * sizeof(int));

    if (axis == 'x') {
        result[0] = x;
        result[1] = (int)(y*cos(a) - z*sin(a));
        result[2] = (int)(y*sin(a) + z*cos(a));
    }
    else if (axis == 'y') {
        result[0] = (int)(x*cos(a) + z*sin(a));
        result[1] = y;
        result[2] = (int)(-x*sin(a) + z*cos(a));
    }
    else if (axis == 'z') {
        result[0] = (int)(x*cos(a) - y*sin(a));
        result[1] = (int)(x*sin(a) + y*cos(a));
        result[2] = z;
    }

    result[0] += pivot[0];
    result[1] += pivot[1];
    result[2] += pivot[2];

    return result;
}

void rotateShape(int** shape, int pointsInShape, int* pivot, char axis, double a) {
    for (int i = 0; i < pointsInShape; i++) {
        int* old_point = shape[i];
        shape[i] = rotatePoint(shape[i], pivot, axis, a);
        free(old_point);
    }
}

int** getCircle(int pointsInCircle, int circle_x0, int circle_y0, int circle_z0, int r) {
    int** circle = (int**)malloc(pointsInCircle * sizeof(int*));
    for (int i = 0; i < pointsInCircle; i++) {
        circle[i] = (int*)malloc(3 * sizeof(int));
    }

    double step = 2*PI/pointsInCircle;
    for (double phi = 0; phi < 2*PI; phi += step) {
        int x = circle_x0 + r * cos(phi);
        int y = circle_y0 + r * sin(phi);

        int pointIndex = phi / step;
        circle[pointIndex][0] = x;
        circle[pointIndex][1] = y;
        circle[pointIndex][2] = circle_z0;
    }

    return circle;
}

void freeShape(int** shape, int pointsInShape) {
    for (int i = 0; i < pointsInShape; i++) {
        free(shape[i]);
    }
    free(shape);
}


int** getDonut(int x0, int y0, int z0, int r0, int r1) {
    int** donut = (int**)malloc(CIRCLES * POINTS_IN_CIRCLE * sizeof(int*));
    for (int i = 0; i < CIRCLES * POINTS_IN_CIRCLE; i++) {
        donut[i] = (int*)malloc(3 * sizeof(int));
    } 

    int** circle = getCircle(POINTS_IN_CIRCLE, x0, y0 - r0, z0, r1);

    int* pivot = (int*)malloc(3 * sizeof(int));
    pivot[0] = x0;
    pivot[1] = y0;
    pivot[2] = z0;

    double step = 2*PI/CIRCLES;
    for (double theta = 0; theta < 2*PI; theta += step) {
        for (int i = 0; i < POINTS_IN_CIRCLE; i++) {
            int index = (int)(theta / step * POINTS_IN_CIRCLE + i);
            int* rotatedPoint = rotatePoint(circle[i], pivot, 'x', theta);
            donut[index][0] = rotatedPoint[0];
            donut[index][1] = rotatedPoint[1];
            donut[index][2] = rotatedPoint[2];
            free(rotatedPoint);
        }
    }
    free(circle);
    free(pivot);
    return donut;
}

void addToTensor(Tensor tensor, int** shape, int pointCount, int x_0, int y_0, int z_0) {
    int pixelArraySize = 12;
    for (int i = 0; i < pointCount; i++) {
        int x = shape[i][0];
        int y = shape[i][1];
        int z = shape[i][2];

        if (x < 0 || x >= SIZE || y < 0 || y >= SIZE || z < 0 || z >= SIZE) {
            continue;
        }

        int vect[3] = {x - x_0, y - y_0, z - z_0};

        double dotProd = vect[0]*lightVect[0] + vect[1]*lightVect[1] + vect[2]*lightVect[2];
        double vectLen = sqrt(vect[0]*vect[0] + vect[1]*vect[1] + vect[2]*vect[2]);
        double lightVectLen = sqrt(lightVect[0]*lightVect[0] + lightVect[1]*lightVect[1] + lightVect[2]*lightVect[2]);
        double angleCos = ((dotProd / vectLen / lightVectLen) + 1) / 2; //0 -> 1

        char pixel = pixels[(int)(round(angleCos * (pixelArraySize-1)))];
        tensor.array[x][y][z] = pixel;
    }
}


Tensor tensor;

int main() {
    tensor = allocateTensor(SIZE, SIZE, SIZE);
    fillTensor(tensor, ' ');
    int** donut = getDonut(SIZE/2, SIZE/2, SIZE/2, SIZE * 0.3, SIZE * 0.15);

    int* pivot = (int*)malloc(3 * sizeof(int));
    pivot[0] = SIZE/2;
    pivot[1] = SIZE/2;
    pivot[2] = SIZE/2;

    srand(time(NULL));

    double angle = PI/2;
    while(1) {
        usleep(50000);
        int** rotatedDonut = (int**)malloc(CIRCLES * POINTS_IN_CIRCLE * sizeof(int*));
        for (int i = 0; i < CIRCLES * POINTS_IN_CIRCLE; i++) {
            rotatedDonut[i] = (int*)malloc(3 * sizeof(int));
            rotatedDonut[i][0] = donut[i][0];
            rotatedDonut[i][1] = donut[i][1];
            rotatedDonut[i][2] = donut[i][2];
        }
        rotateShape(rotatedDonut, CIRCLES * POINTS_IN_CIRCLE, pivot, 'y', angle);
        rotateShape(rotatedDonut, CIRCLES * POINTS_IN_CIRCLE, pivot, 'z', angle);
        fillTensor(tensor, ' ');
        addToTensor(tensor, rotatedDonut, CIRCLES * POINTS_IN_CIRCLE, SIZE/2, SIZE/2, SIZE/2);
        clearScreen();
        printTensorProj(tensor);

        angle += PI/25;
        if (angle > 2*PI)
            angle -= 2*PI;

        freeShape(rotatedDonut, CIRCLES * POINTS_IN_CIRCLE);
    }
    freeTensor(tensor);
    freeShape(donut, CIRCLES * POINTS_IN_CIRCLE);
    return 0;
}

