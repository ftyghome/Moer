/**
 * @file PathIntegrator-new.h
 * @author Chenxi Zhou
 * @brief Path Integrator with new implementation
 * @version 0.1
 * @date 2022-09-21
 * 
 * @copyright NJUMeta (c) 2022 
 * www.njumeta.com
 */

#pragma once

#include "CoreLayer/Ray/Ray.h"
#include "AbstractPathIntegrator.h"
#include "ppg/path.h"
#include "ppg/sdtree.h"

/**
 * @brief Unidirectional path-tracing integrator with new
 * implementation (mitsuba-like)
 * @ingroup Integrator
 */

class PPGIntegrator : public AbstractPathIntegrator 
{
public:
    PPGIntegrator (std::shared_ptr<Camera> _camera,
                       std::unique_ptr<Film> _film,
                       std::unique_ptr<TileGenerator> _tileGenerator,
                       std::shared_ptr<Sampler> _sampler,
                       int _spp,
                       int _renderThreadNum = 4);

    virtual PPGPath getPPGPathByBSDF(const Ray &ray, const std::shared_ptr<Scene> & scene);

    virtual PPGPath getPPGPathBySDTree(const Ray &ray, const std::shared_ptr<Scene> & scene);

    virtual std::pair<Spectrum,PPGPath> LiWithPathByBSDF(const Ray &ray, const std::shared_ptr<Scene> & scene);

    virtual std::pair<Spectrum,PPGPath> LiWithPathBySDTree(const Ray &ray, const std::shared_ptr<Scene> & scene);

    virtual void render(std::shared_ptr<Scene> scene) override;

    void renderPerThread(std::shared_ptr<Scene> scene);

    virtual PathIntegratorLocalRecord evalEmittance(std::shared_ptr<Scene> scene,
                                                    std::optional<Intersection> itsOpt,
                                                    const Ray &ray) override;


    virtual PathIntegratorLocalRecord sampleDirectLighting(std::shared_ptr<Scene> scene,
                                                           const Intersection &its,
                                                           const Ray &ray) override;

    virtual PathIntegratorLocalRecord evalScatter(const Intersection &its,
                                                  const Ray &ray,
                                                  const Vec3d &wi) override;

    virtual PathIntegratorLocalRecord sampleScatter(const Intersection &its,
                                                    const Ray &ray) override;

    virtual double russianRoulette(const Spectrum &T,
                                   int nBounce) override;

    virtual std::pair<std::shared_ptr<Light>, double> chooseOneLight(std::shared_ptr<Scene> scene,
                                                                     double lightSample);

    virtual double chooseOneLightPdf(std::shared_ptr<Scene> scene,
                                     std::shared_ptr<Light> light);

    virtual PathIntegratorLocalRecord evalEnvLights(std::shared_ptr<Scene> scene,
                                                    const Ray &ray);

protected:
    const int nPathLengthLimit = 16;
    const double pRussianRoulette = 0.95;
};