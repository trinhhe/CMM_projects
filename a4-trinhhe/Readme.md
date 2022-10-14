First Name: Henry

Last Name: Trinh

Solution to Question 4:\
[proof](proof.pdf)

Solution to Question 10:\
The optimizer only takes into account the distance between some feature points position to its target position. The handlers' orientation don't matter to the optimization, which is why the handlers always point to right/left side.

Solution to Question 11:\
Our regularizer (reg) has the same size as our control vector u where as the first 4 entries of it are zero and the last 2 entries are used to correct the orientation of our handlers. If right handler moves up, reg[4] has some small positive value to force the right handler to orientate upwards. If right handlers moves down, reg[4] has some small negative number to orientate right handler downwards. Similar case for the left handler.


Assignment writeup: http://crl.ethz.ch/teaching/computational-motion-21/slides/tutorial-a4.pdf

---

Could use ./build.sh on Linux/MacOS
