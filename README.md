# julia-parallel-course-EGU21

#### [vEGU2021: SC4.6 Solving differential equations in parallel with Julia](https://meetingorganizer.copernicus.org/EGU21/session/38986)

👉 **Organisation notes:**
- 💡 All material presented during the short course will be uploaded and available to participants 5 days prior to the event on this repository.
- ⚠️ Live participation to this short course requires [registration to EGU's main meeting](https://egu21.eu/register.html).
- Due to the time-limited course schedule (60min), an interactive overview will cover the course's [objectives](#objectives) instead of the originally planned hands-on.
- ❓ Further interests in solving PDEs in parallel with Julia? Sign-up for a 3h hands-on workshop at [JuliaCon2021](https://juliacon.org/2021/) (free).

----
This short course covers trendy areas in modern geocomputing with broad geoscientific applications. The physical processes governing natural systems' evolution are often mathematically described as systems of differential equations. A performant numerical implementation of the solving algorithm leveraging modern hardware is key and permits to tackle problems that were technically not possible a decade ago.


## Content
* [Objectives](#objectives)


## Objectives
The goal of this short course is to offer an interactive overview on how to solve systems of differential equations in parallel on GPUs using the [Julia language]. [Julia] combines high-level language simplicity and low-level language performance. The resulting codes and applications are fast, short and readable. We will design and implement a numerical algorithm that predicts ice flow dynamics over mountainous topography using a high-performance computing approach.

The online course consists of 2 parts:
1. You will learn about the [Julia language], parallel and distributed computing and iterative solvers.
2. You will implement a PDE solver to predict ice flow dynamics on real topography.

By the end of this short course, you will:
- Have a GPU PDE solver that predicts ice-flow;
- Have a Julia code that achieves similar performance than legacy codes (C, CUDA, MPI);
- Know how the Julia language solves the "two-language problem";
- Be able to leverage the computing power of modern GPU accelerated servers and supercomputers.


[Julia]: https://julialang.org
[Julia language]: https://docs.julialang.org/en/v1/
