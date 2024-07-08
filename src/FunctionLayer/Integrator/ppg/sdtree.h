#pragma once
#include <set>

enum class SpatialSplitType {
    SPLIT_X,
    SPLIT_Y,
    SPLIT_Z
};

struct SDTreeSampleResult {
    Vec3d direction;
    double pdf;
};

inline std::pair<double, double> convert_cart2sphe(const Vec3d &direction) {
    return {fm::acos(direction.z), fm::atan2(direction.y, direction.x)};
}

inline std::pair<double, double> convert_sphe22d(const double theta, const double phi) {
    return {phi / (2 * M_PI), (fm::cos(theta) + 1) / 2};
}

inline std::pair<double, double> convert_cart22d(const Vec3d &direction) {
    const auto [theta, phi] = convert_cart2sphe(direction);
    return convert_sphe22d(theta, phi);
}

inline std::pair<double, double> convert_2d2sphe(const double x, const double y) {
    return {fm::acos(2 * y - 1), x * 2 * M_PI};
}

inline Vec3d convert_sphe2cart(const double theta, const double phi) {
    return {fm::sin(theta) * fm::cos(phi), fm::sin(theta) * fm::sin(phi), fm::cos(theta)};
}

inline Vec3d convert_2d2cart(const double x, const double y) {
    const auto [theta, phi] = convert_2d2sphe(x, y);
    return convert_sphe2cart(theta, phi);
}

struct LRecord {
    Point3d point;
    Spectrum L;
    Vec3d direction;
};

struct DirectionalTreeNode {
    double flux = 0;

    std::shared_ptr<DirectionalTreeNode> quad[4];

    DirectionalTreeNode(const double flux_) : flux(flux_) {}

    DirectionalTreeNode(const DirectionalTreeNode &node) {
        flux = node.flux;
        if (node.quad[0]) {
            for (int i = 0; i < 4; i++) {
                quad[i] = std::make_shared<DirectionalTreeNode>(*node.quad[i]);
            }
        }
    }

    DirectionalTreeNode() = default;

    bool isLeaf() const { return !quad[0]; }
};

struct DirectionalTree {
    std::shared_ptr<DirectionalTreeNode> root;

    DirectionalTree(const std::shared_ptr<DirectionalTree> &orig_dtree, const BoundingBox3f &bbox) {
        root = std::make_shared<DirectionalTreeNode>(*orig_dtree->root);
    }

    DirectionalTree() {
        root = std::make_shared<DirectionalTreeNode>();
    }

    void testSplit(std::shared_ptr<DirectionalTreeNode> &node) {
        assert(!node->quad[0]);
        if (node->flux / root->flux > 0.1) {
            for (auto &i : node->quad) {
                i = std::make_shared<DirectionalTreeNode>(node->flux / 4);
            }
        }
    }

    std::shared_ptr<DirectionalTreeNode> descentToNode(const double x, const double y, const std::function<void(const std::shared_ptr<DirectionalTreeNode> &)> &func) {
        auto cur_node = root;
        func(cur_node);
        Point2d minP{0, 0}, maxP{1, 1};
        if (!cur_node->quad[0]) {
            testSplit(cur_node);
            goto out;
        }
        while (cur_node->quad[0]) {
            Point2d midP = minP + (maxP - minP) * 0.5;
            if (x < midP.x) {
                if (y < midP.y) {
                    cur_node = cur_node->quad[0];
                    coordDescent(minP, maxP, 0);
                } else {
                    cur_node = cur_node->quad[2];
                    coordDescent(minP, maxP, 2);
                }
            } else {
                if (y < midP.y) {
                    cur_node = cur_node->quad[1];
                    coordDescent(minP, maxP, 1);
                } else {
                    cur_node = cur_node->quad[3];
                    coordDescent(minP, maxP, 3);
                }
            }
            assert(cur_node);
            func(cur_node);
            if (!cur_node->quad[0]) {
                testSplit(cur_node);
                break;
            }
        }
    out:
        return cur_node;
    }

    inline static void coordDescent(Point2d &minP, Point2d &maxP, int quadIdx) {
        Point2d midP = minP + (maxP - minP) * 0.5;
        switch (quadIdx) {
            case 0:
                maxP = midP;
                break;
            case 1:
                minP.x = midP.x;
                maxP.y = midP.y;
                break;
            case 2:
                minP.y = midP.y;
                maxP.x = midP.x;
                break;
            case 3:
                minP = midP;
                break;
            default:
                assert(false);
        }
    }

    void apply(const Vec3d &direction_, const double flux) {
        const auto direction = normalize(direction_);
        auto [x, y] = convert_cart22d(direction);

        auto node = descentToNode(x, y, [=](const std::shared_ptr<DirectionalTreeNode> &inode) {
            inode->flux += flux;
        });
    }

    SDTreeSampleResult getSample(double randNum) {
        auto cur_node = root;
        double pdf = 1 / (4 * M_PI);
        Point2d minP{0, 0}, maxP{1, 1};
        while (!cur_node->isLeaf()) {
            int choice = 3;
            double flux[4];
            flux[0] = cur_node->quad[0]->flux;
            for (int i = 1; i < 4; i++) {
                flux[i] = flux[i - 1] + cur_node->quad[i]->flux;
            }
            for (int i = 0; i < 4; i++) {
                if (randNum < flux[i] / flux[3]) {
                    choice = i;
                    break;
                }
            }
            coordDescent(minP, maxP, choice);
            pdf *= 4 * cur_node->quad[choice]->flux / cur_node->flux;
            cur_node = cur_node->quad[choice];
            const double interval = flux[choice] - (choice == 0 ? 0 : flux[choice - 1]);
            randNum = randNum * flux[3] / interval;
            assert(0 <= randNum && randNum <= 1);
        }
        const Point2d sampledP = minP + (maxP - minP) * randNum;
        const Vec3d direction = convert_2d2cart(sampledP.x, sampledP.y);
        return {direction, pdf};
    }
};

struct SpatialTreeNode {
    int numSamples = 0;
    BoundingBox3f bbox{};
    SpatialSplitType split_type = SpatialSplitType::SPLIT_X;
    std::shared_ptr<SpatialTreeNode> l, r;
    std::shared_ptr<DirectionalTree> dtree;
    explicit SpatialTreeNode(const BoundingBox3f &bbox_) : bbox(bbox_) {}
    explicit SpatialTreeNode(const BoundingBox3f &bbox_, SpatialSplitType split_type_, int numSamples_) : numSamples(numSamples_), bbox(bbox_), split_type(split_type_) {}

    std::shared_ptr<SpatialTreeNode> dispatchStep(const Point3d &point) {
        assert(l && r);
        if (l->bbox.contains(point)) return l;
        if (r->bbox.contains(point)) return r;
        return nullptr;
    }

    void dispatchDirectionalTree() {
        assert(l && r);
        l->dtree = std::make_shared<DirectionalTree>(dtree, l->bbox);
        r->dtree = std::make_shared<DirectionalTree>(dtree, r->bbox);
    }

    void split() {
        assert(!l && !r);
        switch (split_type) {
            case SpatialSplitType::SPLIT_X:
                l = std::make_shared<SpatialTreeNode>(BoundingBox3f{bbox.pMin, bbox.pMax - (bbox.pMax - bbox.pMin) * TVector3(.5, .0, .0)}, SpatialSplitType::SPLIT_Y, numSamples / 2);
                r = std::make_shared<SpatialTreeNode>(BoundingBox3f{bbox.pMin + (bbox.pMax - bbox.pMin) * TVector3(.5, .0, .0), bbox.pMax}, SpatialSplitType::SPLIT_Y, numSamples - numSamples / 2);
                break;
            case SpatialSplitType::SPLIT_Y:
                l = std::make_shared<SpatialTreeNode>(BoundingBox3f{bbox.pMin, bbox.pMax - (bbox.pMax - bbox.pMin) * TVector3(.0, .5, .0)}, SpatialSplitType::SPLIT_Z, numSamples / 2);
                r = std::make_shared<SpatialTreeNode>(BoundingBox3f{bbox.pMin + (bbox.pMax - bbox.pMin) * TVector3(.0, .5, .0), bbox.pMax}, SpatialSplitType::SPLIT_Z, numSamples - numSamples / 2);
                break;
            case SpatialSplitType::SPLIT_Z:
                l = std::make_shared<SpatialTreeNode>(BoundingBox3f{bbox.pMin, bbox.pMax - (bbox.pMax - bbox.pMin) * TVector3(.0, .0, .5)}, SpatialSplitType::SPLIT_X, numSamples / 2);
                r = std::make_shared<SpatialTreeNode>(BoundingBox3f{bbox.pMin + (bbox.pMax - bbox.pMin) * TVector3(.0, .0, .5), bbox.pMax}, SpatialSplitType::SPLIT_X, numSamples - numSamples / 2);
                break;
        }
        dispatchDirectionalTree();
        dtree = nullptr;
    }
};

struct SpatialTree {
    int threshold = INT_MAX;
    std::shared_ptr<SpatialTreeNode> root;
    std::shared_ptr<SpatialTreeNode> descentToNode(const Point3d &point, bool incSample = false) {
        auto cur_node = root;
        if (incSample) cur_node->numSamples++;
        while (cur_node->l && cur_node->r) {
            cur_node = cur_node->dispatchStep(point);
            if (incSample) cur_node->numSamples++;
        }
        if (cur_node->numSamples > threshold) {
            cur_node->split();
            cur_node = cur_node->dispatchStep(point);
            assert(cur_node);
        }
        return cur_node;
    }

    void applyVertex(const Vertex &vert_ref, const Vertex &vert_nxt) {
        auto its = vert_ref.itsOpt.value();

        LRecord record;
        record.point = its.position;
        record.direction = vert_nxt.ray.direction;
        record.L = vert_nxt.L * vert_nxt.sampleScatterRecord.f;

        const auto cur_spatial_node = descentToNode(record.point, true);
        assert(cur_spatial_node->dtree);
        const auto cur_directional_node = cur_spatial_node->dtree;
        cur_directional_node->apply(record.direction, record.L.luminance());
    }

    explicit SpatialTree(BoundingBox3f &bbox_) {
        root = std::make_shared<SpatialTreeNode>(bbox_);
        root->dtree = std::make_shared<DirectionalTree>();
    }
};

struct SDTree {
    std::shared_ptr<SpatialTree> stree;
    void applyPath(const PPGPath &path) {
        for (size_t i = 1; i < path.verts.size() - 1; i++) {
            stree->applyVertex(path.verts[i], path.verts[i + 1]);
        }
    }
    void setSpatialSplitThreshold(const int threshold) {
        stree->threshold = threshold;
    }

    SDTreeSampleResult getSampleFromPoint(const Point3d &point, const double randNum) {
        const auto spatial_node = stree->descentToNode(point);
        const auto sample = spatial_node->dtree->getSample(randNum);
        return sample;
    }

    explicit SDTree(BoundingBox3f bbox_) {
        stree = std::make_shared<SpatialTree>(bbox_);
    }
};