# Physics 427 Homework 2

__Due 11:59pm Monday 9/11/2023__

## 1. Root Finding Algorithms

In a C++ header file `problem1.h`, implement the three root-finding algorithms discussed in class. Use template to allow the functions to accept arbitrary functions $f$ for input, and solve the equation $f(x) = 0$. The three functions should look like:

``` c++
template <typename F>
std::tuple<bool, double> fixed_point_iteration(const F& f, double x0, double tolerance);

template <typename F>
std::tuple<bool, double> bisection(const F& f, double a, double b, double tolerance);

template <typename F, typename Fprime>
std::tuple<bool, double> newton_raphson(const F& f, const Fprime& fprime, double x0, double tolerance);
```

Fixed point iteration and Newton's method both start with a trial solution, whereas the bisection method starts with an interval $[a, b]$. All of them accepts a tolerance value. For the purpose of the homework, we take it to mean "how close $f(x)$ is close to 0", meaning that we terminate the iteration when $|f(x_n)| < \epsilon$ where $\epsilon$ is the given tolerance. 

We use `std::tuple` as return value again, since we want to include the information of whether a solution was found. The first element is a boolean value, marking whether the algorithm converged correctly. If the iteration stops due to hitting the tolerance, return `true` for the boolean in the `std::tuple`, and set the second element as the found solution. If for whatever reason the iteration fails to converge, does not hit the target tolerance, or produces NaN, then you should return `false` for the first element in the `std::tuple`, and set the `double` to whatever value the iteration gives you (probably a wrong answer or NaN). Note that in order to achieve this, you may need to set a maximum number of iterations.
    
Note that for Newton's method, you will also need to pass the derivative of $f$ as input, since we need to evaluate it during the iteration.

Finally add and commit your `problem1.h` to the repo. I have included a test C++ file, `test.cpp`, in the repository to help you debug your code. It will assume your functions are implemented in the `problem1.h` header, and try to call these functions to test whether they yield correct results. I also encourage you to write a few of your own test cases, since the homework problem will be graded using a different, more comprehensive set of test cases.

## 2. Finding Lagrange Points

In a system of two celestial bodies orbiting each other, there are equilibrium points for small-mass objects to follow the orbit of the more massive bodies. These are called _Lagrange points_. For example, the recently launched James-Webb Space Telescope (JWST) maintains a small orbit around the second Lagrange point $L_2$ of the Earth-Sun system.

We are interested in the 3 collinear Lagrange points lying on along the line going through the Sun and the Earth. $L_1$ is defined to be the point between the Earth and the Sun where the combined gravitational force of both bodies gives the right centrifugal force for an object to orbit there at the same orbital period as the Earth. This way, such an object will be seen as "stationary" with respect to both the Earth and the Sun.

Let $M_1$ and $M_2$ be the masses of the Sun and the Earth respectively. Let $R$ be the Earth-Sun distance, and $r_1$ be the distance of $L_1$ from the Earth. Both of the Sun and the Earth orbit around their center of mass, which is $RM_2/(M_1 + M_2)$ away from the center of the Sun. Kepler's third law tells us that:

$$
\omega^2 R^3 = G(M_1 + M_2)
$$

where $\omega$ is the orbital frequency, and $G$ is the gravitational constant. An object at $L_1$ will also orbit around the Earth-Sun center of mass, with a radius $r = RM_1/(M_1 + M_2) - r_1$ (see the figure for an illustration). Equating the gravitational force with its centrifugal force at the above orbital frequency, one can write down the equation for $r_1$:

$$
\frac{GM_1}{(R - r_1)^2} - \frac{GM_2}{r_1^2} = \left(\frac{M_1}{M_1 + M_2}R - r_1\right)\frac{M_1 + M_2}{R^3}
$$

Similarly, the second Lagrange point $L_2$ is on the line passing through the Sun and the Earth, except that $L_2$ lies on the far side of the Earth. In this case, an object at $L_2$ orbits at the same frequency as the Earth, but the centrifugal force is given by the combined gravitational force from both the Sun and the Earth. Let $r_2$ be the distance between the center of Earth and $L_2$, a similar reasoning as above gives the equation for $r_2$:

$$
\frac{GM_1}{(R + r_2)^2} + \frac{GM_2}{r_2^2} = \left(\frac{M_1}{M_1 + M_2}R + r_2\right)\frac{M_1 + M_2}{R^3}
$$

The third Lagrange point $L_3$ lies behind the Sun, but the idea is similar to $L_2$. Let the distance between $L_3$ and the center of the Sun be $R - r_3$ (defined this way so that $r_3$ is a small number compared to $R$, just like $r_1$ and $r_2$). The equation for $r_3$ is:

$$
\frac{GM_1}{(R - r_3)^2} + \frac{GM_2}{(2R - r_3)^2} = \left(\frac{M_2}{M_1 + M_2}R + R - r_3\right)\frac{M_1 + M_2}{R^3}
$$

You task is the following: 

First, derive the _dimensionless_ versions of the 3 equations above. You can define $x_{1,2,3} = r_{1,2,3}/R$, and $\mu = M_2/(M_1 + M_2)$ to help with the process. Write down your steps and results in a PDF or Markdown file named `problem2.pdf` or `problem2.md` and commit it to the GitHub repo.

Second, write a C++ program to solve these three _dimensionless_ equations for $x_{1,2,3}$. You should make use of the functions you implemented in Problem 1. Choose the method that makes the most sense to you, and choose a reasonable tolerance. Think what tolerance makes sense: is $10^{-2}$ enough, or do you want to shoot for $10^{-10}$. Print out the answer in the following format:

``` sh
r1: [x1]R
r2: [x2]R
r3: [x3]R
```
For example, if you find $x_1\approx 0.03$, then it should print `r1: 0.03R`. Name your source file `problem2.cpp` and commit it to the GitHub repo.

