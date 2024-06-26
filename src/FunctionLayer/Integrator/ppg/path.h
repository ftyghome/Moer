#pragma once
#include <utility>
#include <vector>

struct Vertex {
    std::optional<Intersection> itsOpt;
    Ray ray;
    PathIntegratorLocalRecord sampleScatterRecord{};
    Spectrum L{.0};
    Vertex(std::optional<Intersection> itsOpt_, const Ray &ray_) : itsOpt(std::move(itsOpt_)), ray(ray_) {}
    Vertex(std::optional<Intersection> itsOpt_, const Ray &ray_, PathIntegratorLocalRecord sampleScatterRecord_)
        : itsOpt(std::move(itsOpt_)), ray(ray_), sampleScatterRecord(std::move(sampleScatterRecord_)) {}
};

struct PPGPath {
    Ray initialRay;
    std::vector<Vertex> verts;
    explicit PPGPath(const Ray &initial_ray) : initialRay(initial_ray) {}
};