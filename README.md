## 🔻ConeCaster

![screenshot](/screenshot.png)

Currently:
* Triangle only
* Single CPU thread
* Regular ray-tracing

Supposed to be:
* Support custom shapes e.g. sdf, volume
* Multi-threaded
* Explore the advantages of cone-tracing

I have tried to implement a BVH for the first time and dived into the topic a bit much. Now I'm learning about BVH generation and traversal methods.

There are 3 different BVH implementations in the [include/Lib/render/scene.hpp](/include/Lib/render/scene.hpp):
* As shown in [Ray Tracing in One Weekend](https://raytracing.github.io/books/RayTracingTheNextWeek.html#boundingvolumehierarchies)
* My interpretation of it into a 8 wide SIMD version
* (🔨WIP) As shown in [Jacco's BVH tutorial series](https://jacco.ompf2.com/2022/04/13/how-to-build-a-bvh-part-1-basics/)

My implementation is the fastest among them, up to ~2x. Though, I haven't finished Jacco's tutorial yet.