#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;

const int WIDTH = 50;
const int HEIGHT = 20;
const double SPEED = 2.0;
const int DELAY_MS = 100;
const int N_SECONDS = 3;

mutex printMutex;

class Rabbit {
public:
    string type;
    double x, y;
    double vx, vy;

    Rabbit(string type, double startX, double startY)
        : type(type), x(startX), y(startY), vx(0), vy(0) {
    }

    void liveAlbino() {
        double currentSpeed = SPEED;

        while (true) {
            {
                lock_guard<mutex> lock(printMutex);
                cout << "[ALBINO] " << type << " at X=" << (int)x << endl;
            }

            x += currentSpeed;

            if (x >= WIDTH || x <= 0) {
                currentSpeed *= -1;
            }

            this_thread::sleep_for(chrono::milliseconds(DELAY_MS));
        }
    }

    void liveOrdinary() {
        changeDirection();

        int stepsForChange = (N_SECONDS * 1000) / DELAY_MS;
        int counter = 0;

        while (true) {
            {
                lock_guard<mutex> lock(printMutex);
                cout << "[ORDINARY] " << type << " at (" << (int)x << ", " << (int)y << ")" << endl;
            }

            x += vx;
            y += vy;

            if (x <= 0 || x >= WIDTH) { x = (x <= 0 ? 0 : WIDTH); vx *= -1; }
            if (y <= 0 || y >= HEIGHT) { y = (y <= 0 ? 0 : HEIGHT); vy *= -1; }

            counter++;
            if (counter >= stepsForChange) {
                {
                    lock_guard<mutex> lock(printMutex);
                    cout << "--- " << type << " changed random direction! ---" << endl;
                }
                changeDirection();
                counter = 0;
            }

            this_thread::sleep_for(chrono::milliseconds(DELAY_MS));
        }
    }

private:
    void changeDirection() {
        double angle = ((double)rand() / RAND_MAX) * 2 * 3.14159;
        vx = SPEED * cos(angle);
        vy = SPEED * sin(angle);
    }
};

int main() {
    srand(time(0));

    vector<thread> threads;

    threads.emplace_back(&Rabbit::liveAlbino, Rabbit("Albino_1", 0, 10));

    threads.emplace_back(&Rabbit::liveOrdinary, Rabbit("Ordinary_A", 25, 10));
    threads.emplace_back(&Rabbit::liveOrdinary, Rabbit("Ordinary_B", 10, 5));

    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }

    return 0;
}