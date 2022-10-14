# Assignment 3 - Boids

## 2 Advanced Time Integration
for h =  1/10\
In explicit Euler particles rotate around the origin faster and the gap between origin and particle widens faster in comparison to the other two methods.\
[explicit Euler with h = 1/10](demos/explicitEuler_1_10.mp4)

In symplectic Euler, we see a few fast full rotations then it slows a bit down. The gaps between origin and particle widens rather slower.\
[symplectic Euler with h = 1/10](demos/symplecticEuler_1_10.mp4)

In explicit Midpoint, rotations may slow down a little bit later than symplectic Euler. Gaps between origin and particles seems kidna similar to symplectic euler.\
[explicitMidpoint with h = 1/10](demos/explicitMidpoint_1_10.mp4)

for h = 1/100\
All 3 methods seem to have similar behaviour.

## Flocking
[cohesion](demos/cohesion.mp4)

[alignment](demos/alignment.mp4)

[separation](demos/separation.mp4)

[collisionAvoidance](demos/collision.mp4)

## Collaborative and Adversarial Behaviors 
[leading](demos/leading.mp4)

The _control\_strategy_ consists of boids staying closer to boids of the same groups (cohesion). Enemy boids radiate a repulsive force (separation). Once a while all boids gather into the mid and try to eat or breed.\
[controlStrategy](demos/control_strategy.mp4)


In addition, I have implemented void _limit\_velocity_ to limit the speed of a boid, _moveToPoint_ moving a boid to a point.

