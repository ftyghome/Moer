/**
 * @file PathIntegrator-new.cpp
 * @author Chenxi Zhou
 * @brief Path Integrator with new implementations
 * @version 0.1
 * @date 2022-09-21
 *
 * @copyright NJUMeta (c) 2022
 * www.njumeta.com
 */

#include "PPGIntergrator.h"
#include "FastMath.h"
#include "ppg/path.h"
#include "FunctionLayer/Material/NullMaterial.h"

#include <thread>

const double eps = 1e-4;

PPGIntegrator::PPGIntegrator(std::shared_ptr<Camera> _camera,
                             std::unique_ptr<Film> _film,
                             std::unique_ptr<TileGenerator> _tileGenerator,
                             std::shared_ptr<Sampler> _sampler,
                             int _spp,
                             int _renderThreadNum) : AbstractPathIntegrator(_camera,
                                                                            std::move(_film),
                                                                            std::move(_tileGenerator),
                                                                            _sampler, _spp,
                                                                            _renderThreadNum) {
}

PPGPath PPGIntegrator::getPPGPathBySDTree(const Ray &initial_ray, const std::shared_ptr<Scene> &scene) {
    return getPPGPathByBSDF(initial_ray,scene);
}


PPGPath PPGIntegrator::getPPGPathByBSDF(const Ray &initial_ray, const std::shared_ptr<Scene> &scene) {

    Ray ray = initial_ray;
    Spectrum throughput{1.0};
    PPGPath ppg_path(ray);
    auto itsOpt = scene->intersect(ray);
    int nBounces = 0;
    PathIntegratorLocalRecord sampleScatterRecord;
    ppg_path.verts.emplace_back(itsOpt, ray, sampleScatterRecord);
    do {
        if (!itsOpt.has_value()) break;
        auto its = itsOpt.value();
        if (its.material->getBxDF(its)->isNull()) {
            ray = Ray{its.position + ray.direction * eps, ray.direction};
            itsOpt = scene->intersect(ray);
            continue;
        }
        nBounces++;

        // Russian roulette.
        double pSurvive = russianRoulette(throughput, nBounces);
        if (sampler->sample1D() >= pSurvive) {
            break;
        }
        // TODO: should I handle this?
        // throughput /= pSurvive;

        sampleScatterRecord = sampleScatter(its, ray);

        if (sampleScatterRecord.f.isBlack() || sampleScatterRecord.pdf == 0) break;

        ray = Ray{its.position + sampleScatterRecord.wi * eps, sampleScatterRecord.wi};
        itsOpt = scene->intersect(ray);
        ppg_path.verts.emplace_back(itsOpt, ray, sampleScatterRecord);
    } while (itsOpt.has_value());
    return ppg_path;
}
std::pair<Spectrum, PPGPath> PPGIntegrator::LiWithPathBySDTree(const Ray &initialRay, const std::shared_ptr<Scene> &scene) {
    Spectrum L{.0};
    auto box = scene->getGlobalBoundingBox();
    Spectrum throughput{1.0};
    auto path = getPPGPathByBSDF(initialRay, scene);
    for (size_t i = 0; i < path.verts.size(); i++) {
        auto &vert = path.verts[i];

        auto [itsOpt, ray, sampleScatterRecord, _] = vert;

        //* ----- The object is emitting light (or just get the env light at infinity) -----
        auto evalLightRecord = evalEmittance(scene, itsOpt, ray);

        if (!evalLightRecord.f.isBlack()) {
            //* The continuous ray hit the emitter or hit nothing but environment lighting.
            //* Multiple importance sampling
            double misw;
            if (i == 0) {
                misw = 1;
            } else {
                misw = MISWeight(sampleScatterRecord.pdf, evalLightRecord.pdf);
                if (sampleScatterRecord.isDelta) {
                    //* MIS will not be applied with delta distribution of BSDF.
                    misw = 1.0;
                }
            }

            vert.L = evalLightRecord.f * misw;
        }
        if (!itsOpt.has_value()) break;
        auto its = itsOpt.value();

        //* ----- Direct lighting to the object, and reflect to the direction -----
        for (int i = 0; i < nDirectLightSamples; ++i) {
            PathIntegratorLocalRecord sampleLightRecord = sampleDirectLighting(scene, its, ray);
            PathIntegratorLocalRecord evalScatterRecord = evalScatter(its, ray, sampleLightRecord.wi);

            if (!sampleLightRecord.f.isBlack()) {
                //* Multiple importance sampling
                double misw = MISWeight(sampleLightRecord.pdf, evalScatterRecord.pdf);
                if (sampleLightRecord.isDelta) {
                    // * MIS will not be applied with delta distribution of light source.
                    misw = 1.0;
                }
                vert.L += sampleLightRecord.f * evalScatterRecord.f / sampleLightRecord.pdf * misw / nDirectLightSamples;
            }
        }
    }

    // Backward light collection

    for (size_t i = path.verts.size() - 1; i >= 1; i--) {
        auto &cur = path.verts[i], &prev = path.verts[i - 1];
        prev.L += cur.L * cur.sampleScatterRecord.f / cur.sampleScatterRecord.pdf;
    }

    return {path.verts.empty() ? 0 : path.verts[0].L, path};
}

std::pair<Spectrum, PPGPath> PPGIntegrator::LiWithPathByBSDF(const Ray &initialRay,
                                                             const std::shared_ptr<Scene> &scene) {
    Spectrum L{.0};
    auto box = scene->getGlobalBoundingBox();
    Spectrum throughput{1.0};
    auto path = getPPGPathByBSDF(initialRay, scene);
    for (size_t i = 0; i < path.verts.size(); i++) {
        auto &vert = path.verts[i];

        auto [itsOpt, ray, sampleScatterRecord, _] = vert;

        //* ----- The object is emitting light (or just get the env light at infinity) -----
        auto evalLightRecord = evalEmittance(scene, itsOpt, ray);

        if (!evalLightRecord.f.isBlack()) {
            //* The continuous ray hit the emitter or hit nothing but environment lighting.
            //* Multiple importance sampling
            double misw;
            if (i == 0) {
                misw = 1;
            } else {
                misw = MISWeight(sampleScatterRecord.pdf, evalLightRecord.pdf);
                if (sampleScatterRecord.isDelta) {
                    //* MIS will not be applied with delta distribution of BSDF.
                    misw = 1.0;
                }
            }

            vert.L = evalLightRecord.f * misw;
        }
        if (!itsOpt.has_value()) break;
        auto its = itsOpt.value();

        //* ----- Direct lighting to the object, and reflect to the direction -----
        for (int i = 0; i < nDirectLightSamples; ++i) {
            PathIntegratorLocalRecord sampleLightRecord = sampleDirectLighting(scene, its, ray);
            PathIntegratorLocalRecord evalScatterRecord = evalScatter(its, ray, sampleLightRecord.wi);

            if (!sampleLightRecord.f.isBlack()) {
                //* Multiple importance sampling
                double misw = MISWeight(sampleLightRecord.pdf, evalScatterRecord.pdf);
                if (sampleLightRecord.isDelta) {
                    // * MIS will not be applied with delta distribution of light source.
                    misw = 1.0;
                }
                vert.L += sampleLightRecord.f * evalScatterRecord.f / sampleLightRecord.pdf * misw / nDirectLightSamples;
            }
        }
    }

    // Backward light collection

    for (size_t i = path.verts.size() - 1; i >= 1; i--) {
        auto &cur = path.verts[i], &prev = path.verts[i - 1];
        prev.L += cur.L * cur.sampleScatterRecord.f / cur.sampleScatterRecord.pdf;
    }

    return {path.verts.empty() ? 0 : path.verts[0].L, path};
}
void PPGIntegrator::render(std::shared_ptr<Scene> scene) {
    std::vector<std::thread> threads;
    for (int i = 0; i < renderThreadNum; i++) {
        threads.emplace_back(&PPGIntegrator::renderPerThread, this, scene);
    }

    for (int i = 0; i < renderThreadNum; i++) {
        threads[i].join();
    }

    printProgress(1.f);
}
void PPGIntegrator::renderPerThread(std::shared_ptr<Scene> scene) {
    /**
     * @warning Other part of Integrator uses the original sampler
     *          so I have to fill its vectors here.
     *          In fact every time a fresh sampler is needed we should
     *          use Sampler::clone() to get one.
     */
    static int tileFinished = 0;

    sampler->startPixel({0, 0});
    auto ssampler = sampler->clone(0);
    while (true) {
        auto optionalTile = tileGenerator->generateNextTile();
        if (optionalTile == std::nullopt)
            break;
        auto tile = optionalTile.value();

        for (auto it = tile->begin(); it != tile->end(); ++it) {

            auto bbox_ = scene->getGlobalBoundingBox();
            BoundingBox3f bbox{Point3d{
                                   bbox_.pMin.x - .1,
                                   bbox_.pMin.y - .1,
                                   bbox_.pMin.z - .1,
                               },
                               Point3d{
                                   bbox_.pMax.x + .1,
                                   bbox_.pMax.y + .1,
                                   bbox_.pMax.z + .1}};

            SDTree cur_tree(bbox), nxt_tree(bbox);

            auto pixelPosition = *it;

            const auto &cam = *this->camera;
            /**
             * @bug Sampler is NOT designed for multi-threads, need copy for each thread.
             *      Sampler::clone() will return a Sampler copy, only with same sampling
             *      strategy, random numbers are not guaranteed to be identical.
             */
            // sampler->startPixel(pixelPosition);
            ssampler->startPixel(pixelPosition);

            cur_tree.setSpatialSplitThreshold(12);

            for (int i = 0; i < 1; i++) {
                auto [L, path] = LiWithPathByBSDF(
                    cam.generateRay(
                        film->getResolution(),
                        pixelPosition,
                        ssampler->getCameraSample()),
                    scene);
                // std::cout << "Applying path " << std::endl;
                cur_tree.applyPath(path);
                // std::cout << "Applying path complete" << std::endl;
            }
            for (int i = 0; i < spp; i++) {
                auto L = LiWithPathBySDTree(
                    cam.generateRay(
                        film->getResolution(),
                        pixelPosition,
                        ssampler->getCameraSample()),
                    scene).first;
                film->deposit(pixelPosition, L);
                /**
                 * @warning spp used in this for loop belongs to Integrator.
                 *          It is irrelevant with spp passed to Sampler.
                 *          And Sampler has no sanity check for subscript
                 *          of sample vector. Error may occur if spp passed to
                 *          Integrator is bigger than which passed to Sampler.
                 */
                ssampler->nextSample();
            }
        }

        //* Finish one tile rendering
        if (++tileFinished % 5) {
            printProgress((float)tileFinished / tileGenerator->tileCount);
        }
    }
}

/// @brief Eval surface or infinite light source (on itsOpt) radiance and ignore medium transmittance.
/// @param scene Scene description. Used to query scene lighting condition.
/// @param itsOpt Current intersection point which could be a light source. If there's no intersection, eval the radiance of environment light.
/// @param ray Current ray which connects last intersection point and itsOpt.
/// @return Current ray direction, obtained light radiance and solid angle dependent pdf. Note that there is no corresponding sampling process for pdf and the pdf value should NOT be applied to calculate the final radiance contribution.
PathIntegratorLocalRecord PPGIntegrator::evalEmittance(std::shared_ptr<Scene> scene,
                                                       std::optional<Intersection> itsOpt,
                                                       const Ray &ray) {
    Vec3d wo = -ray.direction;
    Spectrum LEmission(0.0);
    double pdfDirect = 1.0;
    if (!itsOpt.has_value()) {
        auto record = evalEnvLights(scene, ray);
        LEmission = record.f;
        pdfDirect = record.pdf;
    } else if (itsOpt.value().object && itsOpt.value().object->getLight()) {
        auto its = itsOpt.value();
        Normal3d n = its.geometryNormal;
        auto light = itsOpt.value().object->getLight();
        auto record = light->eval(ray, its, ray.direction);
        LEmission = record.s;
        Intersection tmpIts;
        tmpIts.position = ray.origin;
        pdfDirect = record.pdfDirect * chooseOneLightPdf(scene, light);
    }
    // Path integrator will ignore medium transmittance. And (of course) there will be no occlusion.
    // Spectrum transmittance(1.0);
    return {ray.direction, LEmission, pdfDirect, false};
}

/// @brief Sample on the distribution of direct lighting and ignore medium transmittance.
/// @param scene Scene description. Multiple shadow ray intersect operations will be performed.
/// @param its Current intersection point which returned pdf dependent on.
/// @param ray Current ray. Should only be applied for time records.
/// @return Sampled direction on the distribution of direct lighting and corresponding solid angle dependent pdf. An extra flag indicites that whether it sampled on a delta distribution.
PathIntegratorLocalRecord PPGIntegrator::sampleDirectLighting(std::shared_ptr<Scene> scene,
                                                              const Intersection &its,
                                                              const Ray &ray) {
    auto [light, pdfChooseLight] = chooseOneLight(scene, sampler->sample1D());
    auto record = light->sampleDirect(its, sampler->sample2D(), ray.timeMin);
    double pdfDirect = record.pdfDirect * pdfChooseLight;// pdfScatter with respect to solid angle
    Vec3d dirScatter = record.wi;
    Spectrum Li = record.s;
    Point3d posL = record.dst;
    Point3d posS = its.position;
    Spectrum transmittance(1.0);// todo: transmittance eval
    Ray visibilityTestingRay(posL - dirScatter * 1e-4, -dirScatter, ray.timeMin, ray.timeMax);
    auto visibilityTestingIts = scene->intersect(visibilityTestingRay);
    if (!visibilityTestingIts.has_value() || visibilityTestingIts->object != its.object || (visibilityTestingIts->position - posS).length2() > 1e-6) {
        transmittance = 0.0;
    }
    if (!visibilityTestingIts.has_value() && light->lightType == ELightType::INFINITE) {
        transmittance = 1.0;
    }
    return {dirScatter, Li * transmittance, pdfDirect, record.isDeltaPos};
}

/// @brief Eval the scattering function, i.e., bsdf * cos. The cosine term will not be counted for delta distributed bsdf (inside f).
/// @param its Current intersection point which is used to obtain local coordinate and bxdf.
/// @param ray Current ray which is used to calculate bsdf $f(\omega_i,\omega_o)$.
/// @param dirScatter (already sampled) scattering direction.
/// @return scattering direction, bsdf value and bsdf pdf.
PathIntegratorLocalRecord PPGIntegrator::evalScatter(const Intersection &its,
                                                     const Ray &ray,
                                                     const Vec3d &dirScatter) {
    if (its.material != nullptr) {
        std::shared_ptr<BxDF> bxdf = its.material->getBxDF(its);
        Normal3d n = its.geometryNormal;
        double wiDotN = fm::abs(dot(n, dirScatter));
        Vec3d wi = its.toLocal(dirScatter);
        Vec3d wo = its.toLocal(-ray.direction);
        return {
            dirScatter,
            bxdf->f(wo, wi, false) * wiDotN,
            bxdf->pdf(wo, wi),
            false};
    } else {
        // Path integrator will igonre phase function distribution.
        return {};
    }
}

/// @brief Sample a direction along with a pdf value in term of solid angle according to the distribution of bsdf.
/// @param its Current intersection point.
/// @param ray Current incident ray.
/// @return Sampled scattering direction, bsdf * cos, corresponding pdf and whether it is sampled on a delta distribution.
PathIntegratorLocalRecord PPGIntegrator::sampleScatter(const Intersection &its,
                                                       const Ray &ray) {
    if (its.material != nullptr) {
        Vec3d wo = its.toLocal(-ray.direction);
        std::shared_ptr<BxDF> bxdf = its.material->getBxDF(its);
        Vec3d n = its.geometryNormal;
        BxDFSampleResult bsdfSample = bxdf->sample(wo, sampler->sample2D(), false);
        double pdf = bsdfSample.pdf;
        Vec3d dirScatter = its.toWorld(bsdfSample.directionIn);
        double wiDotN = fm::abs(dot(dirScatter, n));
        return {dirScatter, bsdfSample.s * wiDotN, pdf, BxDF::MatchFlags(bsdfSample.bxdfSampleType, BXDF_SPECULAR)};
    } else {
        // todo: sample phase function
        return {};
    }
}

/// @brief Russian roulette method.
/// @param throughput Current thorughput, i.e., multiplicative (bsdf * cos) / pdf.
/// @param nBounce Current bounce depth.
/// @return Survive probility after Russian roulette.
double PPGIntegrator::russianRoulette(
    const Spectrum &throughput,
    int nBounce) {
    // double pSurvive = std::min(pRussianRoulette, throughput.sum());
    double pSurvive = pRussianRoulette;
    if (nBounce > nPathLengthLimit)
        pSurvive = 0.0;
    if (nBounce <= 20)
        pSurvive = 1.0;
    return pSurvive;
}

/// @brief (discretely) sample a light source.
/// @param scene Scene description which is used to query scene lights.
/// @param lightSample A random number within [0,1].
/// @return Pointer of the sampled light source and corresponding (discrete) probility.
std::pair<std::shared_ptr<Light>, double>
PPGIntegrator::chooseOneLight(std::shared_ptr<Scene> scene,
                              double lightSample) {
    // uniformly weighted
    std::shared_ptr<std::vector<std::shared_ptr<Light>>> lights = scene->getLights();
    int numLights = lights->size();
    int lightID = std::min(numLights - 1, (int)(lightSample * numLights));
    std::shared_ptr<Light> light = lights->operator[](lightID);
    return {light, 1.0 / numLights};
}

/// @brief Calculate the (discrete) probility that a specific light source is sampled.
/// @param scene Scene description which is used to query scene lights.
/// @param light The specific light source.
/// @return Corresponding (discrete) probility.
double PPGIntegrator::chooseOneLightPdf(std::shared_ptr<Scene> scene,
                                        std::shared_ptr<Light> light) {
    std::shared_ptr<std::vector<std::shared_ptr<Light>>> lights = scene->getLights();
    int numLights = lights->size();
    return 1.0 / numLights;
}

/// @brief Eval the effect of environment light source (infinite area light).
/// @param scene Scene description which is used to query scene lights.
/// @param ray Current ray.
/// @return light source direction, light radiance and pdf. Without sampling process, the pdf can NOT be used to calculate final radiance.
PathIntegratorLocalRecord PPGIntegrator::evalEnvLights(std::shared_ptr<Scene> scene,
                                                       const Ray &ray) {
    std::shared_ptr<std::vector<std::shared_ptr<Light>>> lights = scene->getLights();
    Spectrum L(0.0);
    double pdf = 0.0;
    for (auto light : *lights) {
        auto record = light->evalEnvironment(ray);
        L += record.s;
        pdf += record.pdfEmitDir;
    }
    return {-ray.direction, L, pdf};
}
