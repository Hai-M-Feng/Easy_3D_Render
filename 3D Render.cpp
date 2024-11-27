#include <cstdio>
#include <iostream>
#include <graphics.h>
#include <thread>
#include <conio.h>
#include <vector>
#include <random>
#include <ctime>
#include <algorithm>
#include <cmath>

using namespace std;

const double M_PI = 3.14159265358979; //PI
const double ScreenWidth = 500; //屏幕宽
const double ScreenHeight = 500; //屏幕高

///////////////////////////////////////////////////////////////////////////////////////  点相关

//三维点或向量类
class Point_3D {
public:
    double p_X;
    double p_Y;
    double p_Z;
    double p_W;

    // 构造方法
    Point_3D(double x, double y, double z, double w) : p_X(x), p_Y(y), p_Z(z), p_W(w) {}
    Point_3D() : Point_3D(0, 0, 0, 0) {}

    bool inScreen() {
        return ((p_X >= 0 && p_X <= ScreenWidth) && (p_Y >= 0 && p_Y <= ScreenHeight));
    }

    // 重构向量和小数的乘法
    Point_3D operator*(const double& d) const {
        return Point_3D(p_X * d, p_Y * d, p_Z * d, p_W * d);
    }

    // 重构向量和向量的加法
    Point_3D operator+(const Point_3D& other) const {
        double newX = this->p_X + other.p_X;
        double newY = this->p_Y + other.p_Y;
        double newZ = this->p_Z + other.p_Z;
        double newW = this->p_W + other.p_W;

        // 规范化结果，确保 w 为 1
        if (newW != 0) {
            return Point_3D(newX / newW, newY / newW, newZ / newW, 1.0);
        }
        else {
            return Point_3D(newX, newY, newZ, 1.0); // 或者返回其他合理的默认值
        }
    }

    // 重构向量和向量的减法
    Point_3D operator-(const Point_3D& other) const {
        double newX = this->p_X - other.p_X;
        double newY = this->p_Y - other.p_Y;
        double newZ = this->p_Z - other.p_Z;
        double newW = this->p_W - other.p_W;

        // 规范化结果，确保 w 为 1
        if (newW != 0) {
            return Point_3D(newX / newW, newY / newW, newZ / newW, 1.0);
        }
        else {
            return Point_3D(newX, newY, newZ, 1.0); // 或者返回其他合理的默认值
        }
    }
};

// 计算向量的模
double magnitude(const Point_3D& v) {
    return sqrt(v.p_X * v.p_X + v.p_Y * v.p_Y + v.p_Z * v.p_Z);
}

//用来归一化向量
Point_3D normalize(const Point_3D& vec) {
    double len = magnitude(vec);
    if (len != 0) {
        return Point_3D(vec.p_X / len, vec.p_Y / len, vec.p_Z / len, vec.p_W);
    }
    else {
        return vec;
    }
}

//计算两个向量的叉积
Point_3D crossProduct(const Point_3D& v1, const Point_3D& v2) {
    return Point_3D(
        v1.p_Y * v2.p_Z - v1.p_Z * v2.p_Y,
        v1.p_Z * v2.p_X - v1.p_X * v2.p_Z,
        v1.p_X * v2.p_Y - v1.p_Y * v2.p_X,
        0.0
    );
}

///////////////////////////////////////////////////////////////////////////////////////  矩阵相关

//矩阵类，四阶
class Matrix {
public:
    vector<vector<double>> matrix;
    Matrix() : matrix(4, vector<double>(4, 0)) {}

    //向量和矩阵的乘法
    Point_3D operator*(const Point_3D& point) const {
        double x = matrix[0][0] * point.p_X + matrix[0][1] * point.p_Y + matrix[0][2] * point.p_Z + matrix[0][3] * point.p_W;
        double y = matrix[1][0] * point.p_X + matrix[1][1] * point.p_Y + matrix[1][2] * point.p_Z + matrix[1][3] * point.p_W;
        double z = matrix[2][0] * point.p_X + matrix[2][1] * point.p_Y + matrix[2][2] * point.p_Z + matrix[2][3] * point.p_W;
        double w = matrix[3][0] * point.p_X + matrix[3][1] * point.p_Y + matrix[3][2] * point.p_Z + matrix[3][3] * point.p_W;

        // 归一化结果，确保 w 为 1
        if (w != 0) {
            return Point_3D(x / w, y / w, z / w, 1.0);
        }
        else {
            return Point_3D(x, y, z, 1.0); // 或者返回其他合理的默认值
        }
    }

    //矩阵之间的乘法
    Matrix operator*(const Matrix& other) const {
        Matrix result;
        for (int i = 0; i < 4; ++i) {  //遍历结果矩阵的行
            for (int j = 0; j < 4; ++j) {  //遍历结果矩阵的列
                for (int k = 0; k < 4; ++k) {  //遍历中间维度
                    result.matrix[i][j] += matrix[i][k] * other.matrix[k][j];
                }
            }
        }
        return result;
    }
};

//创建旋转矩阵，需要输入摄像机的方向向量
Matrix createRotationMatrix(const Point_3D& direction) {
    // 确保方向向量已经规范化
    Point_3D zAxis = normalize(direction);

    // 选择一个初始向量作为 x 轴
    Point_3D xAxis;
    if (fabs(zAxis.p_X) < 0.99 && fabs(zAxis.p_Y) < 0.99) {
        xAxis = normalize(Point_3D(1, 0, 0, 0));
    }
    else {
        xAxis = normalize(Point_3D(0, 1, 0, 0));
    }

    // 计算 y 轴
    Point_3D yAxis = normalize(crossProduct(zAxis, xAxis));

    // 重新计算 x 轴，确保正交性
    xAxis = normalize(crossProduct(yAxis, zAxis));

    // 构建旋转矩阵
    Matrix rotationMatrix;
    rotationMatrix.matrix = {
        {xAxis.p_X, yAxis.p_X, zAxis.p_X, 0},
        {xAxis.p_Y, yAxis.p_Y, zAxis.p_Y, 0},
        {xAxis.p_Z, yAxis.p_Z, zAxis.p_Z, 0},
        {0, 0, 0, 1}
    };

    return rotationMatrix;
}

//创建相机平移矩阵，应该输入相机的世界坐标
Matrix createTranslationalMatrix(const Point_3D& position) {
    Matrix translationalMatrix;
    translationalMatrix.matrix = {
        {1, 0, 0, -position.p_X},
        {0, 1, 0, -position.p_Y},
        {0, 0, 1, -position.p_Z},
        {0, 0, 0, 1}
    };

    return translationalMatrix;
}

//创建相机的方向矩阵，需要输入相机的位置和方向向量
Matrix createOrientationMatrix(const Point_3D& position, const Point_3D& direction) {
    Point_3D zAxis = normalize(direction);
    Point_3D upVector = (fabs(zAxis.p_X) < 0.99 && fabs(zAxis.p_Y) < 0.99) ? normalize(Point_3D(0, 1, 0, 0)) : normalize(Point_3D(1, 0, 0, 0));
    Point_3D rightVector = normalize(crossProduct(zAxis, upVector));
    upVector = normalize(crossProduct(rightVector, zAxis));

    Matrix orientationMatrix;
    orientationMatrix.matrix = {
        {rightVector.p_X, upVector.p_X, zAxis.p_X, 0},
        {rightVector.p_Y, upVector.p_Y, zAxis.p_Y, 0},
        {rightVector.p_Z, upVector.p_Z, zAxis.p_Z, 0},
        {0, 0, 0, 1}
    };

    return orientationMatrix;
}

// 创建透视投影矩阵
Matrix createPerspectiveMatrix(double fov, double aspectRatio, double nearPlane, double farPlane) {
    double f = 1.0 / tan((fov * M_PI / 180.0) / 2.0);
    Matrix perspectiveMatrix;
    perspectiveMatrix.matrix = {
        {f / aspectRatio, 0, 0, 0},
        {0, f, 0, 0},
        {0, 0, (farPlane + nearPlane) / (nearPlane - farPlane), (2 * farPlane * nearPlane) / (nearPlane - farPlane)},
        {0, 0, -1, 0}
    };

    return perspectiveMatrix;
}

//创建裁剪矩阵
Matrix createClipMatrix(double fov, double aspectRatio, double nearPlane, double farPlane) {
    double f = 1.0 / tan((fov * M_PI / 180.0) / 2.0);
    Matrix clipMatrix;
    clipMatrix.matrix = {
        {f / aspectRatio, 0, 0, 0},
        {0, f, 0, 0},
        {0, 0, (farPlane + nearPlane) / (nearPlane - farPlane), (2 * farPlane * nearPlane) / (nearPlane - farPlane)},
        {0, 0, -1, 0}
    };

    return clipMatrix;
}

//创建屏幕矩阵
Matrix createScreenMatrix(double screenWidth, double screenHeight) {
    Matrix screenMatrix;
    screenMatrix.matrix = {
        {(screenWidth / 2), 0, 0, (screenWidth / 2)},
        {0, -(screenHeight / 2), 0, (screenHeight / 2)},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };

    return screenMatrix;
}

///////////////////////////////////////////////////////////////////////////////////////  线相关

//线类
class Line {
public:
    Point_3D p1, p2; //线的两个点
    bool ISrend; //是否渲染

    //构造方法
    Line(Point_3D p1t, Point_3D p2t) : p1(p1t), p2(p2t) {
        ISrend = (p1.inScreen() || p2.inScreen());
    }

    //setter
    void setLine(Point_3D& p1t, Point_3D& p2t) {
        p1 = p1t;
        p2 = p2t;
        ISrend = (p1.inScreen() || p2.inScreen());
    }

    //返回线的两个点
    vector<Point_3D> getPoints() {
        vector<Point_3D> points;
        points.push_back(p1);
        points.push_back(p2);
        return points;
    }
};

void drawLine(Line& l) {
    if (l.ISrend) {
        line(l.p1.p_X, l.p1.p_Y, l.p2.p_X, l.p2.p_Y);
        printf("line drawed at %.2lf,%.2lf~%.2lf,%.2lf\n", l.p1.p_X, l.p1.p_Y, l.p2.p_X, l.p2.p_Y);
    }
    else {
        printf("line out of screen!\n");
    }
}

///////////////////////////////////////////////////////////////////////////////////////

//正方体类
class Cube {
public:
    unsigned Hash_ID; //用于区分每个正方体 
    Point_3D center; //正方体的中心位置
    Point_3D direction; //正方体的朝向
    double length; //正方体的边长

    //构造方法
    Cube(Point_3D p, Point_3D d, double l) : center(p), direction(d), length(l) { Hash_ID = rand(); }

    //返回正方体的所有边
    vector<Line> getLines() {
        // 首先计算正方体的半边长
        double halfLength = length / 2;

        // 计算正方体的8个顶点
        vector<Point_3D> vertices(8);

        // 创建旋转矩阵
        Matrix rotationMatrix = createRotationMatrix(direction);

        // 计算原始未旋转状态下的顶点位置
        Point_3D baseVertices[] = {
            Point_3D(-halfLength, -halfLength, -halfLength, 1),
            Point_3D(halfLength, -halfLength, -halfLength, 1),
            Point_3D(halfLength, halfLength, -halfLength, 1),
            Point_3D(-halfLength, halfLength, -halfLength, 1),
            Point_3D(-halfLength, -halfLength, halfLength, 1),
            Point_3D(halfLength, -halfLength, halfLength, 1),
            Point_3D(halfLength, halfLength, halfLength, 1),
            Point_3D(-halfLength, halfLength, halfLength, 1)
        };

        // 应用旋转和平移变换到每个顶点
        for (int i = 0; i < 8; ++i) {
            Point_3D rotated = rotationMatrix * baseVertices[i];
            vertices[i] = Point_3D(rotated.p_X + center.p_X, rotated.p_Y + center.p_Y, rotated.p_Z + center.p_Z, 1);
        }

        // 定义正方体的12条棱
        vector<Line> lines;
        int edges[][2] = {
            {0, 1}, {1, 2}, {2, 3}, {3, 0},  // 下底面
            {4, 5}, {5, 6}, {6, 7}, {7, 4},  // 上底面
            {0, 4}, {1, 5}, {2, 6}, {3, 7}   // 连接上下底面的边
        };

        // 根据顶点创建线
        for (auto& edge : edges) {
            lines.emplace_back(vertices[edge[0]], vertices[edge[1]]);
        }

        return lines;
    }
};

//摄像机类
class Camera {
public:
    Point_3D position; //摄像机的位置
    Point_3D direction; //摄像机的朝向
    double FOV = 60.0; //摄像机的视野
    double screenRate = ScreenWidth / ScreenHeight;

    //构造方法
    Camera(Point_3D p, Point_3D d) : position(p), direction(d) {};
        
    //移动相机
    void move(double x, double y, double z) {
        position.p_X += x;
        position.p_Y += y;
        position.p_Z += z;
    }

    // 绕 X 轴旋转方向向量（世界坐标系）
    void rotateX(double angle) {
        double radian = angle * M_PI / 180.0;
        Matrix rotationMatrix;
        rotationMatrix.matrix = {
            {1, 0, 0, 0},
            {0, cos(radian), -sin(radian), 0},
            {0, sin(radian), cos(radian), 0},
            {0, 0, 0, 1}
        };
        direction = rotationMatrix * direction;
    }

    // 绕 Y 轴旋转方向向量（世界坐标系）
    void rotateY(double angle) {
        double radian = angle * M_PI / 180.0;
        Matrix rotationMatrix;
        rotationMatrix.matrix = {
            {cos(radian), 0, sin(radian), 0},
            {0, 1, 0, 0},
            {-sin(radian), 0, cos(radian), 0},
            {0, 0, 0, 1}
        };
        direction = rotationMatrix * direction;
    }

    // 绕 Z 轴旋转方向向量（世界坐标系）
    void rotateZ(double angle) {
        double radian = angle * M_PI / 180.0;
        Matrix rotationMatrix;
        rotationMatrix.matrix = {
            {cos(radian), -sin(radian), 0, 0},
            {sin(radian), cos(radian), 0, 0},
            {0, 0, 1, 0},
            {0, 0, 0, 1}
        };
        direction = rotationMatrix * direction;
    }

    // 移动视角
    void moveView(double dx, double dy, double dz) {
        if (dx != 0) {
            rotateX(dx);
        }
        if (dy != 0) {
            rotateY(dy);
        }
        if (dz != 0) {
            rotateZ(dz);
        }
    }

    //渲染画面 测试
    void render(vector<Cube> list) {
        //创建矩阵
        Matrix T = createTranslationalMatrix(position);
        Matrix O = createOrientationMatrix(position, direction);
        Matrix P = createPerspectiveMatrix(FOV, screenRate, 1.0, 1000.0);
        Matrix C = createClipMatrix(FOV, screenRate, 1.0, 1000.0);
        Matrix S = createScreenMatrix(ScreenWidth, ScreenHeight);
        Matrix FinalMatrix = S * C * P * O * T;

        vector<Line> lineListBefore;
        vector<Line> lineListAfter;

        //导出正方形所有线
        for (Cube& tempCube : list) {
            vector<Line> tempLines = tempCube.getLines();
            lineListBefore.insert(lineListBefore.end(), tempLines.begin(), tempLines.end());
        }

        // 处理所有线
        for (Line& tempLine : lineListBefore) {
            vector<Point_3D> tempPoints = tempLine.getPoints();
            Point_3D a = FinalMatrix * tempPoints[0];
            Point_3D b = FinalMatrix * tempPoints[1];

            printf("transform line %.2lf,%.2lf~%.2lf,%.2lf ", tempLine.p1.p_X, tempLine.p1.p_Y, tempLine.p2.p_X, tempLine.p2.p_Y);
            printf("to %.2lf,%.2lf~%.2lf,%.2lf\n", a.p_X, a.p_Y, b.p_X, b.p_Y);
            lineListAfter.push_back(Line(a, b));
        }

        //绘制所有线
        BeginBatchDraw();
        cleardevice();
        for (Line& tempLine : lineListAfter) {
            drawLine(tempLine);
        }
        EndBatchDraw();
    }
};

int main() {
    srand(unsigned(time(NULL) * 100));

    vector<Cube> render_list;
    render_list.push_back(Cube({ 0, 0, 0, 1 }, { 0, 0, -1, 0 }, 100));
    //render_list.push_back(Cube({ 0, 0, 150, 1 }, { 0, 0, -1, 0 }, 50));
    Camera camera_1({ 0, 0, 200 ,1 }, { 0, 0, -1, 0 });

    initgraph(ScreenWidth, ScreenHeight);

    camera_1.render(render_list);

    while (true) {
        char c = _getch();
        switch (c) {
        case 'w':
        case 'W': {
            // 向前移动（沿世界Z轴）
            camera_1.move(0, 0, -10);
            camera_1.render(render_list);
            break;
        }
        case 's':
        case 'S': {
            // 向后移动（沿世界Z轴）
            camera_1.move(0, 0, 10);
            camera_1.render(render_list);
            break;
        }
        case 'a':
        case 'A': {
            // 向左移动（沿世界X轴）
            camera_1.move(-10, 0, 0);
            camera_1.render(render_list);
            break;
        }
        case 'd':
        case 'D': {
            // 向右移动（沿世界X轴）
            camera_1.move(10, 0, 0);
            camera_1.render(render_list);
            break;
        }
        case 'f':
        case 'F': {
            // 向顶端移动（沿世界y轴）
            camera_1.move(0, 10, 0);
            camera_1.render(render_list);
            break;
        }
        case 'c':
        case 'C': {
            // 向底端移动（沿世界y轴）
            camera_1.move(0, -10, 0);
            camera_1.render(render_list);
            break;
        }
        case 'h':
            camera_1.moveView(0, 10, 0); //向左转摄像机视角
            camera_1.render(render_list);
            break;
        case 'k':
            camera_1.moveView(0, -10, 0); //向右转摄像机视角
            camera_1.render(render_list);
            break;
        case 'u':
            camera_1.moveView(10, 0, 0); //向上转摄像机视角
            camera_1.render(render_list);
            break;
        case 'j':
            camera_1.moveView(-10, 0, 0); //向下转摄像机视角
            camera_1.render(render_list);
            break;
        }
    }



    closegraph();
    
    return 0;
}