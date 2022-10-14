#pragma once

#include <sim/RB.h>
#include <sim/RBRenderer.h>
#include <sim/RBSpring.h>
#include <sim/RBUtils.h>

using namespace std;
namespace crl {

// subfunctions
Quaternion updateRotationGivenAngularVelocity(const Quaternion &q,
                                              const V3D &angularVelocity,
                                              double dt) {
    double angularVelocityMagnitude = angularVelocity.norm();
    // avoid divide by close to zero...
    if (angularVelocityMagnitude > 1e-10) {
        // TODO: Ex.1 Integration
        // implement quaternion update logic
        // q_p = rot(w, dt) * q
        Quaternion rot;
        rot.w() = cos(dt*angularVelocityMagnitude/2.0);
        rot.vec() = angularVelocity/angularVelocityMagnitude * sin(dt*angularVelocityMagnitude/2.0);
        // TODO: fix the following line.
        return rot * q;
    }
    return q;
}

/**
 * our simulation world governed by the rigid-body dynamics
 */
class RBPhysicsEngine {
public:
    // constructor
    RBPhysicsEngine() {}

    // desctructor
    ~RBPhysicsEngine() {
        // we are going to delete every rigid body when simulation world is
        // destroyed.
        for (uint i = 0; i < rbs.size(); i++) delete rbs[i];
        rbs.clear();
        for (uint i = 0; i < springs.size(); i++) delete springs[i];
        springs.clear();
    }

    /**
     * add rigid body to simulation world.
     * note that we assume this is a cube block with approx. 0.24 x 0.24 x 0.24. 
     */
    RB *addRigidBodyToEngine() {
        rbs.push_back(new RB());
        // we use default mass = 100 kg
        rbs.back()->rbProps.mass = 100;
        rbs.back()->rbProps.collision = false;
        rbs.back()->rbProps.id = rbs.size() - 1;
        // add force and torque
        f_ext.push_back(V3D());
        tau_ext.push_back(V3D());
        return rbs.back();
    }

    /**
     * add rigid body with collision to simulation world.
     * note that we assume this is a sphere with radius = 0.1. 
     */
    RB *addCollidingRigidBodyToEngine() {
        double i = 0.4 * 100 * 0.1 * 0.1;
        rbs.push_back(new RB());
        // we use default mass = 100 kg
        rbs.back()->rbProps.mass = 100;
        rbs.back()->rbProps.setMOI(i, i, i, 0, 0, 0);
        rbs.back()->rbProps.collision = true;
        rbs.back()->rbProps.id = rbs.size() - 1;
        // add force and torque
        f_ext.push_back(V3D());
        tau_ext.push_back(V3D());
        return rbs.back();
    }

    /**
     * add spring to simulation world.
     */
    RBSpring *addSpringToEngine(RB *parent, RB *child, P3D pJPos, P3D cJPos) {
        springs.push_back(new RBSpring());
        springs.back()->parent = parent;
        springs.back()->child = child;
        // local position of attached point from parent/child frames
        springs.back()->pJPos = pJPos;
        springs.back()->cJPos = cJPos;
        // we use default spring constant = 2000;
        springs.back()->k = 10000;

        // default rest length is the distance between attaching points when
        // the spring is added.
        if (parent == nullptr) {
            // TODO: Ex.2-1
            // implement your logic for a spring which is attached to world
            //
            // TODO: Fix the following line
            springs.back()->l0 = V3D(pJPos - child->state.getWorldCoordinates(cJPos)).norm();
        } else {
            // TODO: Ex.2-2
            // implement your logic for a spring where both ends are attached
            // to rigid bodies
            //
            // TODO: Fix the following line
            springs.back()->l0 = V3D(parent->state.getWorldCoordinates(pJPos) - child->state.getWorldCoordinates(cJPos)).norm();
        }
        return springs.back();
    }

    /**
     * apply external force (no spring force, no gravity. Force comes from 
     * third sources) to rigid body.
     */
    void applyForceTo(RB *rb, const V3D &f, const P3D &p) {
        // add force only if rb is in rbs
        for (uint i = 0; i < rbs.size(); i++) {
            if (rbs[i] == rb) {
                f_ext[i] += f;
                V3D r = rb->state.getWorldCoordinates(V3D(p));
                tau_ext[i] += r.cross(f);
            }
        }
    }

    /**
     * simulation stepping logic. advance one simulation timestep with dt.
     * the unit of dt is second. 
     */
    void step(double dt) {
        // external force and torque
        for (uint i = 0; i < rbs.size(); i++) {
            // force and torque by gravity
            f_ext[i] += rbs[i]->rbProps.mass * V3D(0, RBGlobals::g, 0);
        }

        // force and torque by springs
        for (RBSpring *spring : springs) {
            // TODO: Ex.2 Spring force
            // compute spring force f_spring = -kx and torque tau_spring and
            // add them to f_ext and tau_ext
            //
            // Hint:
            // - spring->l0 is the rest length and spring->k is the spring contant
            // - you can retrieve index of rb in this->rbs list from rb->rbProps.id

            if (spring->parent == nullptr) {
                // TODO: Ex.2-1
                // implement your logic for a spring which is attached to world
                P3D childpin_world = spring->child->state.getWorldCoordinates(spring->cJPos);
                V3D vect = V3D(childpin_world, spring->pJPos);
                double current_length = vect.norm();
                double displacement = spring->l0 - current_length;
                // std::cout << "displacement: " << displacement << std::endl;
                V3D spring_force((-1*spring->k * displacement)*vect.normalized());
                f_ext[spring->child->rbProps.id] += spring_force;
                tau_ext[spring->child->rbProps.id] += V3D(spring->child->state.pos, childpin_world).cross(spring_force);
            } else {
                // TODO: Ex.2-2
                // implement your logic for a spring where both ends are attached
                // to rigid bodies.
                P3D childpin_world = spring->child->state.getWorldCoordinates(spring->cJPos);
                P3D parentpin_world = spring->parent->state.getWorldCoordinates(spring->pJPos);
                V3D vect = V3D(parentpin_world - childpin_world);
                double current_length = vect.norm();
                double displacement = spring->l0 - current_length;
                V3D spring_force((-1*spring->k * displacement)*vect.normalized());
                f_ext[spring->child->rbProps.id] += spring_force;
                tau_ext[spring->child->rbProps.id] += V3D(spring->child->state.pos, childpin_world).cross(spring_force);
                f_ext[spring->parent->rbProps.id] -= spring_force;
                tau_ext[spring->parent->rbProps.id] -= V3D(spring->parent->state.pos, parentpin_world).cross(spring_force);
            }
        }

        for (uint i = 0; i < rbs.size(); i++) {
            RB *rb = rbs[i];
            // retrieve saved force and tau
            V3D f = f_ext[i];
            V3D tau = tau_ext[i];

            // TODO: Ex.1 Integration
            // implement forward (explicit) Euler integration scheme for computing velocity and pose.
            //
            // Hint:
            // - recall, you need to compute 3x3 moment of inertia matrix expressed in world frame
            //
            // implement your logic here.
            // double mass = rb->rbProps.mass;
            // Matrix3x3 rot = rb->state.orientation.toRotationMatrix();
            // Matrix3x3 MOI_world = rot * rb->rbProps.MOI_local * rot.transpose();
            // V3D vel_i = rb->state.velocity;
            // V3D ang_i = rb->state.angularVelocity;
            // rb->state.velocity = vel_i + dt*f/mass;
            // rb->state.angularVelocity = ang_i + dt*MOI_world.inverse()*(tau - ang_i.cross((Vector3d) (MOI_world * ang_i)));
            // rb->state.pos = rb->state.pos + dt*vel_i;
            // rb->state.orientation = updateRotationGivenAngularVelocity(rb->state.orientation, ang_i, dt);

            // TODO: Ex.3 Stable Simulation
            // why our simulation is blown up? let's make it more stable!
            //
            // comment out (do not erase!) your logic for Ex.1 and implement Ex.3 here.

            //symplectic euler
            double mass = rb->rbProps.mass;
            Matrix3x3 rot = rb->state.orientation.toRotationMatrix();
            Matrix3x3 MOI_world = rot * rb->rbProps.MOI_local * rot.transpose();
            rb->state.velocity = rb->state.velocity + dt*f/mass;
            rb->state.angularVelocity = rb->state.angularVelocity + dt*MOI_world.inverse()*(tau - rb->state.angularVelocity.cross((Vector3d) (MOI_world * rb->state.angularVelocity)));
            //pos and orientation update after next if statement
           
            if (simulateCollisions && rb->rbProps.collision) {
                // TODO: Ex.4 Impulse-based Collisions
                // we will simulate collisions between a spherical rigidbody and
                // the ground plane. implement impulse-based collisions here. use
                // coefficient of restituation "epsilon". (it's a member variable 
                // of this class). we assume the friction is infinite.
                // note that we ignore collisions between rigidbodies.
                //
                // Steps:
                // 0. read the material "ImpulseBasedCollisions" on CMM21 website
                // carefully.
                // 1. update linear and angular velocity (not pose yet!!) of the 
                // rigidbody by external force and torque before this if statement 
                // block
                // 2. if rb->rbProps.collision == true, we will assume this rb is
                // a spherical object, and collide with the ground. (if not, 
                // it's a cube and it does not collide with any other objects.)
                // 3. compute impulse
                // 4. update linear and angular velocity with impulse
                // 5. now update pose of the rigid body (by integrating velocity)
                //
                // Hint:
                // - the radius of the sphere is 0.1 m
                // - detect collision if 1) the y coordinate of the point at the
                // bottom of the sphere < 0 and 2) the y component of linear
                // velocity of the point at the botton < 0.
                // - we will assume that a collision only happens at the bottom
                // points.
                // - we will assume there's only one contact between a sphere
                // and the ground

                //
                // Ex.4 implementation here
                //

                //check if collision
                P3D bottom_point = rb->state.pos - P3D(0,0.1,0);
                if(bottom_point[1] < 0 && rb->state.getVelocityForPoint_global(bottom_point)[1] < 0) {
                    V3D r = V3D(rb->state.pos, bottom_point);
                    V3D u_rel = rb->state.velocity + rb->state.angularVelocity.cross(r);
                    V3D contact_normal = V3D(bottom_point, rb->state.pos).normalized();

                    Matrix3x3 K = Matrix3x3::Identity()/mass - getCrossProductMatrix(r) * MOI_world * getCrossProductMatrix(r);

                    V3D J(K.inverse()*(-u_rel - eps*u_rel[1]*contact_normal));
                    cout << "J: " << J << endl;
                    cout << "mass " << mass << endl;
                    rb->state.velocity = rb->state.velocity + J/mass;
                    rb->state.angularVelocity = rb->state.angularVelocity + MOI_world.inverse()*(r.cross(J));
                }
                

            }

            rb->state.pos = rb->state.pos + dt*rb->state.velocity;
            rb->state.orientation = updateRotationGivenAngularVelocity(rb->state.orientation, rb->state.angularVelocity, dt);
        }

        // clean up
        for (uint i = 0; i < f_ext.size(); i++) {
            f_ext[i] = V3D();
            tau_ext[i] = V3D();
        }
    }

    /**
     * draw every rigid body belongs to world.
     */
    inline void draw(const gui::Shader &rbShader) {
        // draw moi boxes
        for (uint i = 0; i < this->rbs.size(); i++) {
            if (!this->rbs[i]->rbProps.fixed) {
                if (this->rbs[i]->rbProps.collision)
                    crl::RBRenderer::drawCollisionRB(this->rbs[i], rbShader);
                else
                    crl::RBRenderer::drawMOI(this->rbs[i], rbShader);
            }
        }

        // draw springs
        for (uint i = 0; i < this->springs.size(); i++) {
            P3D start, end;
            if (this->springs[i]->parent == nullptr) {
                start = this->springs[i]->pJPos;
                end = this->springs[i]->child->state.getWorldCoordinates(
                    this->springs[i]->cJPos);
            } else {
                start = this->springs[i]->parent->state.getWorldCoordinates(
                    this->springs[i]->pJPos);
                end = this->springs[i]->child->state.getWorldCoordinates(
                    this->springs[i]->cJPos);
            }
            drawCylinder(start, end, 0.05, rbShader);
        }

        // and now coordinate frames
        if (showCoordFrame) {
            for (uint i = 0; i < this->rbs.size(); i++)
                crl::RBRenderer::drawCoordFrame(this->rbs[i], rbShader);
        }
    }

    /**
     * returns NULL if no RBs are hit by the ray...
     */
    RB *getFirstRBHitByRay(const Ray &ray, P3D &intersectionPoint) {
        RB *selectedRB = nullptr;
        double t = DBL_MAX;
        P3D tmpIntersectionPoint = P3D(0, 0, 0);

        for (uint i = 0; i < rbs.size(); i++) {
            if (rbs[i]->getRayIntersectionPoint(ray, tmpIntersectionPoint)) {
                double tTmp = ray.getRayParameterFor(tmpIntersectionPoint);
                if (tTmp < t) {
                    selectedRB = rbs[i];
                    t = tTmp;
                    intersectionPoint = tmpIntersectionPoint;
                }
            }
        }
        return selectedRB;
    }

public:
    // this is a list of all rigid bodies and springs belong to the world.
    std::vector<RB *> rbs;
    std::vector<RBSpring *> springs;

    // coefficients
    float eps = 0.0;  // restitution

    // drawing flags
    bool showCoordFrame = true;

    // options
    bool simulateCollisions = false;

private:
    // list of force and tau applied to rigid bodies
    std::vector<V3D> f_ext;
    std::vector<V3D> tau_ext;
};
}  // namespace crl