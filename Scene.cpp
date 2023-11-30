#include "Scene.h"

Scene::Scene() {}
Scene::Scene(const Scene& scene) {}
Scene::~Scene() {}

Scene* Scene::snowman() {
	Scene* scene = new Scene();

	Shape* snowball = Shape::generateSnowball(Vector2f(1, .16), .15, Vector2f(0, 0));
	scene->snow_shapes.push_back(snowball);

	Shape* snowball2 = Shape::generateSnowball(Vector2f(1, .4), .1, Vector2f(0, 0));
	scene->snow_shapes.push_back(snowball2);

	Shape* snowball3 = Shape::generateSnowball(Vector2f(1, .57), .07, Vector2f(0, 0));
	scene->snow_shapes.push_back(snowball3);

	Shape* snowball4 = Shape::generateSnowball(Vector2f(1.8, .4), .08, Vector2f(-500, -40));
	scene->snow_shapes.push_back(snowball4);

	return scene;
}

Scene* Scene::snowballSmash() {

	Scene* scene = new Scene();

	Shape* snowball = Shape::generateSnowball(Vector2f(1.6, .6), .15, Vector2f(-200, 0));
	scene->snow_shapes.push_back(snowball);

	Shape* snowball2 = Shape::generateSnowball(Vector2f(.3, .45), .15, Vector2f(200, 0));
	scene->snow_shapes.push_back(snowball2);

	return scene;
}

Scene* Scene::highspeedSnowballSmash() {
	Scene* scene = new Scene();

	Shape* snowball = Shape::generateSnowball(Vector2f(1.6, .5), .15, Vector2f(-400, 0));
	scene->snow_shapes.push_back(snowball);

	Shape* snowball2 = Shape::generateSnowball(Vector2f(.3, .5), .15, Vector2f(400, 0));
	scene->snow_shapes.push_back(snowball2);

	return scene;
}
