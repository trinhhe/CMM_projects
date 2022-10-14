#ifndef BOIDS_H
#define BOIDS_H
#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/Sparse>
#include <functional>
#include <cmath>
template <typename T, int dim>
using Vector = Eigen::Matrix<T, dim, 1, 0, dim, 1>;
using namespace std;

template <typename T, int n, int m>
using Matrix = Eigen::Matrix<T, n, m, 0, n, m>;

// add more for yours
enum MethodTypes {
        FREEFALL=0, SEPARATION=1, ALIGNMENT=2, COHESION=3, LEADER=4, COLLISIONOBJ=5, 
        COLLABORATIVE=6
    };

template <class T, int dim>
class Boids
{
    typedef Matrix<T, Eigen::Dynamic, 1> VectorXT;
    typedef Matrix<T, dim,Eigen::Dynamic> TVStack;
    typedef Vector<T, dim> TV;
    typedef Matrix<T, dim, dim> TM;
    
private:
    TVStack positions;
    TVStack velocities;
    int n;
    bool update = false;
    T dt = 1.0/80;
    T M;
    TV mouse_position;
    //low/high limits for random velocity generation
    T LO = -0.3;
    T HI = 0.3;
    T max_speed = 1;
    // boids belonging to grp 0 or 1
    Eigen::VectorXf grp;
    int breeding_time;
    int eating_time;
    int war_time;

public:
    Boids() :n(1) {}
    Boids(int n) :n(n) {
        initializePositions();
    }
    ~Boids() {}

    void setParticleNumber(int n) {n = n;}
    int getParticleNumber() { return n; }
    void initializePositions()
    {
        // cout << "INIT POS" << endl;
        positions = TVStack::Zero(dim, n).unaryExpr([&](T dummy){return static_cast <T> (rand()) / static_cast <T> (RAND_MAX);}); 
        // positions.col(0) = TV(0,0);
        // positions.col(1) = TV(1, 1);
        velocities = TVStack::Zero(dim, n).unaryExpr([&](T dummy){return LO + static_cast <T> (rand()) /static_cast <T> (RAND_MAX/(HI-LO));}); 
        M = 1;
        breeding_time = 0;
        eating_time = 0;
        war_time = 0;
        grp = Eigen::VectorXf::Zero(n);
        grp.tail(n/2) = Eigen::VectorXf::Constant(n/2, 1);
        // cout << "pos: " << endl;
        // cout << positions << " ";
        // cout << endl << "------------------------------------" << endl;

        // initiliaze velocities for circular motion
        // for(int i = 0; i < velocities.cols(); i++){
        //     velocities(0, i) = -positions(1,i);
        //     velocities(1, i) = positions(0,i);
        // }
    }

    void setMousePosition(TV x)
    {
        mouse_position = x;
    }

    T getVelocity(TV v)
    {
        T temp = 0;
        for(int i = 0; i < dim; i++) {
            temp += v(i) * v(i);
        }
        return sqrt(temp);
    }

    // limit magnitude of velocity to avoid overshoot
    void limit_velocity()
    {
        for(int i = 0; i < velocities.cols(); i++) {
            T vel = getVelocity(velocities.col(i));
            if (vel > max_speed)
                velocities.col(i) = (velocities.col(i) / vel) * max_speed;
        }
    }

    void updateBehavior(MethodTypes type)
    {
        if(!update)
            return;

        switch (type)
        {  
            case FREEFALL:
                symplecticEuler([this](TVStack pos){return gravity(pos);});
                break;

            case SEPARATION:
                symplecticEuler([this](TVStack pos){return separation(pos);});
                break;

            case ALIGNMENT:
                symplecticEuler([this](TVStack pos){return alignment(pos);});
                break;

            case COHESION:
                symplecticEuler([this](TVStack pos){return cohesion(pos);});
                break;

            case LEADER:
                symplecticEuler([this](TVStack pos){return flocking_leading(pos);});
                break;

            case COLLISIONOBJ:
                symplecticEuler([this](TVStack pos){return flocking_leading_avoiding(pos);});
                // symplecticEuler([this](TVStack pos){return noforce(pos);});

                break;
            
            case COLLABORATIVE:
                symplecticEuler([this](TVStack pos){return control_strategy(pos);});
                if(breeding_time % 500 == 0)
                    breeding();
                breeding_time++;
                if(eating_time % 500 == 0)
                    eating();
                eating_time++;

                // cout << positions.cols() << endl;
                // cout << "breed time" <<breeding_time << endl;
                break;

            // case EXPLIEULER:
            //     explicitEuler([this](TVStack pos){return circle(pos);});
            //     break;
            // case SYMPLIEULER:
            //     symplecticEuler([this](TVStack pos){return circle(pos);});
            //     break;
            // case EXPLIMIDPOINT:
            //     explicitMidpoint([this](TVStack pos){return circle(pos);});
            //     break;
            default:
                break;
        }
    }
    void pause()
    {
        update = !update;
    }
    TVStack getPositions()
    {
        return positions;
    }

    VectorXT getGrp() 
    {
        return grp;
    }

    T getVector(TV p)
    {
        T temp = 0;
        for(int i = 0; i < dim; i++) {
            temp += p(i) * p(i);
        }
        return sqrt(temp);
    }
    
    //distance between p and q
    T getDistance(TV p, TV q)
    {
        T temp = 0;
        for(int i = 0; i < dim; i++) {
            temp += (p(i)-q(i)) * (p(i)-q(i));
        }
        return sqrt(temp);
    }

    bool checkNeighbourhood(TV p, TV center, T radius)
    {
        T temp = 0;
        for(int i = 0; i < dim; i++) {
            temp += (p(i)-center(i))*(p(i)-center(i));
        }
        return (temp <= radius*radius);        
    }

    void explicitEuler(std::function<TVStack (TVStack)> f) 
    {
        TVStack current_velocities = velocities;
        velocities = current_velocities + dt * M * f(positions);
        positions = positions + dt * current_velocities;
    }

    void symplecticEuler(std::function<TVStack (TVStack)> f)
    {
        positions = positions + dt * velocities;
        velocities = velocities + dt * M * f(positions);
        // cout << "velocities before limit: " << endl;
        // cout << velocities << endl;
        limit_velocity();
        // cout << "velocities after limit: " << endl;
        // cout << velocities << endl;
    }

    void explicitMidpoint(std::function<TVStack (TVStack)> f)
    {
        TVStack positions_half = positions + dt/2.0 * velocities;
        TVStack velocities_half = velocities + dt/2.0 * M * f(positions);
        positions = positions + dt * velocities_half;
        velocities = velocities + dt * M * f(positions_half);
    }

    TVStack noforce(TVStack positions) 
    {
        return TVStack::Zero(dim, n);
    }

    TVStack gravity(TVStack positions) 
    {
        TVStack force(positions.rows(), positions.cols());
        for(int i = 0; i < positions.cols(); i++)
            force.col(i) = Eigen::Vector2f(0,0.2);
        return force;
    }

    // circle motion around origin
    TVStack circle(TVStack positions)
    {
        TVStack force(positions.rows(), positions.cols());
        for(int i = 0; i < positions.cols(); i++){
            T radius = getVector(positions.col(i));
            T particle_velocity = sqrt(velocities(0, i) * velocities(0,i) + velocities(1,i) * velocities(1,i));
            T angular_vel = particle_velocity/radius;
            T x = -positions(0, i) * angular_vel * angular_vel;
            T y = -positions(1, i) * angular_vel * angular_vel;
            force.col(i) = Eigen::Vector2f(x,y);
        }
        return force;
    }

    TVStack cohesion(TVStack positions)
    {
        TVStack force(positions.rows(), positions.cols());
        force.setZero();
        T cohension_range = 0.4;
        T numberNeighbours = 0;
        TV averagePosition;
        averagePosition.setZero();
        for(int i = 0; i < positions.cols(); i++) {
            for(int j = 0; j < positions.cols(); j++) {
                //only consider other birds
                if(j != i){
                    if(checkNeighbourhood(positions.col(j), positions.col(i), cohension_range)) {
                        numberNeighbours++;
                        averagePosition += positions.col(j);
                    }
                }
            }
            if(numberNeighbours != 0) {
                averagePosition /= numberNeighbours;
                force.col(i) = (averagePosition - positions.col(i));
                averagePosition.setZero();
                numberNeighbours = 0;
            }
        }
        return force;
    }

    TVStack alignment(TVStack positions)
    {
        TVStack force(velocities.rows(), velocities.cols());
        force.setZero();
        T alignment_range = 0.1;
        T numberNeighbours = 0;
        TV averageVelocity;
        averageVelocity.setZero();
        for(int i = 0; i < velocities.cols(); i++) {
            for(int j = 0; j < velocities.cols(); j++) {
                //only consider other birds
                if(j != i){
                    if(checkNeighbourhood(positions.col(j), positions.col(i), alignment_range)) {
                        numberNeighbours++;
                        averageVelocity += velocities.col(j);
                    }
                }
            }
            if(numberNeighbours != 0) {
                averageVelocity /= numberNeighbours;
                // force.col(i) = (averageVelocity - velocities.col(i)).normalized();
                // force.col(i) = (averageVelocity - velocities.col(i));
                force.col(i) = averageVelocity - velocities.col(i);
                averageVelocity.setZero();
                numberNeighbours = 0;
            }
        }
        return force;
    }

    TVStack separation(TVStack positions)
    {
      TVStack force(positions.rows(), positions.cols());
        force.setZero();
        // T numberNeighbours = 0;
        T separation_range = 0.1;
        TV offset;
        offset.setZero();
        for(int i = 0; i < positions.cols(); i++) {
            for(int j = 0; j < positions.cols(); j++) {
                //only consider other birds
                if(j != i){
                    if(checkNeighbourhood(positions.col(j), positions.col(i), separation_range)) {
                        T dist = getDistance(positions.col(i), positions.col(j));
                        //if they're too close add some manual offset to it
                        if (dist < 0.0125) {
                            offset += (positions.col(i) - positions.col(j)) + TV(0.1,0.1);
                        } else {
                            offset += (positions.col(i) - positions.col(j));
                        }
                        // numberNeighbours++;
                    }
                }
            }
            
            force.col(i) = offset;
            offset.setZero();
                // numberNeighbours = 0;
            
        }
        return force;   
    }

    TVStack moveToPoint(TVStack positions, TV point)
    {
        TVStack force(positions.rows(), positions.cols());
        force.setZero();
        for(int i = 0; i < velocities.cols(); i++) {
            force.col(i) = point - positions.col(i);
        }
        return force;   
    }
    
    TVStack flocking_leading(TVStack positions)
    {
        TVStack force(positions.rows(), positions.cols());
        force.setZero();
        TVStack f1 = cohesion(positions);
        TVStack f2 = alignment(positions);
        TVStack f3 = leading(positions);
        TVStack f4 = separation(positions);

        //damping coefficents
        T a = 0.3;
        T b = 0.4;
        T c = 3;
        T d = 15;
        
        force = a*f1 + b*f2 + c*f3 + d*f4;
        // force.col(0) = c*f3.col(0);
        return force;
    }

     TVStack flocking_leading_avoiding(TVStack positions)
    {
        TVStack force(positions.rows(), positions.cols());
        force.setZero();
        TVStack f1 = cohesion(positions);
        TVStack f2 = alignment(positions);
        TVStack f3 = leading(positions);
        TVStack f4 = separation(positions);
        TVStack f5 = avoidCollision(positions);

        //damping coefficents
        T a = 0.3;
        T b = 0.4;
        T c = 5;
        T d = 7;
        T e = 10;
        force = a*f1 + b*f2 + c*f3 + d*f4 + e*f5;
        // force.col(0) = c*f3.col(0);
        return force;
    }

    TVStack avoidCollision(TVStack positions)
    {
        TVStack force(positions.rows(), positions.cols());
        force.setZero();
        TV circle_center = TV(0.5, 0.5);
        T circle_radius = 0.17;
        for(int i = 0; i < velocities.cols(); i++){
            //position where boid is headed to, multiplied by a factor to increase boids view range
            TV ahead = (positions.col(i) + velocities.col(i))*50;
            //check if collsion with circle
            T distance = getDistance(ahead, circle_center);
            //check if boid is in circle
            T distance1 = getDistance(positions.col(i), circle_center);  
            // boid is too close to object, add repulsive force directly to velocity vector
            if(distance < circle_radius || distance1 < circle_radius){
                if(positions.col(i)[0] < 0.5 && positions.col(i)[1] < 0.5)
                    force.col(i) = -(ahead - circle_center);
                else if (positions.col(i)[0] >= 0.5 && positions.col(i)[1] >= 0.5)
                    force.col(i) = (ahead - circle_center);
                else if (positions.col(i)[0] >= 0.5 && positions.col(i)[1] < 0.5) {
                    force.col(i) = (ahead - circle_center);
                    force.col(i)[1] *= -1;
                }
                else {
                    force.col(i) = (ahead - circle_center);
                    force.col(i)[0] *= -1;
                }

            }
        }
        return force;
    }

    TVStack leading(TVStack positions)
    {
        TVStack force(positions.rows(), positions.cols());
        force.setZero();
        //force for leader following mouse position
        force.col(0) = mouse_position - positions.col(0);
        //rest are following leader
        for(int i = 1; i < positions.cols(); i++) {
            force.col(i) = positions.col(0) - positions.col(i);
        }
        return force;
    }

    void breeding()
    {
        T breeding_radius = 0.05;
        T currentNumberBoids = positions.cols();
        // openToBreeding(i) = 0 means that broid i is still open to breed in this time iteration
        VectorXT openToBreeding = VectorXT::Zero(currentNumberBoids);
        for(int i = 0; i < currentNumberBoids; i++) {
            for(int j = 0; j < currentNumberBoids; j++) {
                if(i!=j && grp(i) == grp(j) && openToBreeding(i) == 0 && openToBreeding(j) == 0) {
                    T distance = getDistance(positions.col(i), positions.col(j));
                    if(distance < breeding_radius) {
                        TVStack pos(dim, 1);
                        pos.col(0) = TV::Zero(dim).unaryExpr([&](T dummy){return static_cast <T> (rand()) / static_cast <T> (RAND_MAX);});
                        TVStack vel(dim, 1);
                        vel.col(0) =  TV::Zero(dim).unaryExpr([&](T dummy){return LO + static_cast <T> (rand()) /static_cast <T> (RAND_MAX/(HI-LO));}); 
                        n++;
                        positions.conservativeResize(dim, n);
                        velocities.conservativeResize(dim, n);
                        grp.conservativeResize(n, 1);
                        positions.col(n-1) = pos;
                        velocities.col(n-1) = vel;
                        grp(n-1) = grp(j);
                        openToBreeding(i) = 1;
                        openToBreeding(j) = 1;
                    }
                }
            }
        }
    }

    void eating()
    {
        T eating_radius = 0.05;
        // hungry(i) = 0 means that boid i is still hungry
        VectorXT hungry = VectorXT::Zero(positions.cols());
        //hold index of 3 three boids
        Eigen::Vector3i index(0,0,0);
        int count = 0;
        int currentNumberBoids = positions.cols() - 1;
        for(int i = 0; i < positions.cols(); i++) {
            for(int j = 0; j < positions.cols(); j++) {
                if(i!=j && grp(i) != grp(j) && hungry(j) == 0 && count < 3) {
                    T distance = getDistance(positions.col(i), positions.col(j));
                    if(distance < eating_radius) {
                        index(count) = j;
                        count++;
                    }
                }
            }
            if(count == 3){
                //boid i to be eaten
                positions.col(i) = positions.col(currentNumberBoids);
                velocities.col(i) = velocities.col(currentNumberBoids);
                grp(i) = grp(currentNumberBoids);
                hungry(index(0)) = 1;
                hungry(index(1)) = 1;
                hungry(index(2)) = 1;
                count = 0;
                currentNumberBoids--;
            }
        }
        n = currentNumberBoids + 1;
        positions.conservativeResize(dim, currentNumberBoids + 1);
        velocities.conservativeResize(dim, currentNumberBoids + 1);
        grp.conservativeResize(currentNumberBoids + 1, 1);
        // cout << "after positions: " << positions.cols() << endl;

    }


    TVStack control_strategy (TVStack dummy)
    {
        TVStack force(positions.rows(), positions.cols());
        force.setZero();
        TVStack f1 = cohesion(positions);
        TVStack f2 = alignment(positions);
        TVStack f3 = separation(positions);
        // TVStack f4 = moveToPoint(positions, TV(0.5,0.5));
        double danger_radius = 0.5;
        double stay_radius = 0.3;
        T average = 0;
        T average1 = 0;
        TV cummu = TV(0,0);
        TV cummu1 = TV(0,0);
        for(int i = 0; i < positions.cols(); i++) {
            for(int j = 0; j < positions.cols(); j++) {
                if(j != i){
                    //cohesion for same grp
                    T dist = getDistance(positions.col(i), positions.col(j));
                    if(dist < stay_radius && grp(i) == grp(j)){
                        cummu += positions.col(j) - positions.col(i);
                        // average++;
                    }
                    //separation for different grp
                    if(dist < danger_radius && grp(i) != grp(j)){
                        cummu1 -= (positions.col(j) - positions.col(i));
                        // average1++;
                    }
                }
            }
            // force.col(i) = cummu/average;
            // force.col(i) += cummu1/average1;
            force.col(i) = cummu;
            force.col(i) += cummu1;
        }
        war_time++;
        for(int i = 0; i < positions.cols(); i++) {
            if(war_time % 1000 == 0) {
                force.col(i) = 50* (TV(0.5,0.5) - positions.col(i));
            } else {
                force.col(i) += TV(0,0) - positions.col(i);            
                force.col(i) += TV::Zero(dim).unaryExpr([&](T dummy){return -10 + static_cast <T> (rand()) /static_cast <T> (RAND_MAX/(10+10));});
            }      
        }
        
        //damping coefficents
        T a = 3;
        T b = 5.5;
        T c = 7;
        force +=  a*f1 + b*f2 + c*f3;
        return force;
    }


};
#endif
