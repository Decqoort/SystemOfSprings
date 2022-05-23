#include <iostream>
#include <SFML/Graphics.hpp>
#include <fstream>
#include <iomanip>

class Body {
public:
	double x, v, a;
	Body() {
		x = 0;
		v = 0;
		a = 0;
	}
	Body(double x, double v) :x(x), v(v), a(0){}
	Body operator+(const Body& a)const {
		return Body{ this->x + a.x, this->v + a.v };
	}
	Body operator-(const Body& a)const {
		return Body{ this->x - a.x, this->v - a.v };
	}
	Body operator*(const Body& a)const {
		return Body{ this->x * a.x, this->v * a.v };
	}
	Body operator*(const double a)const {
		return Body{ this->x * a, this->v * a };
	}
	Body operator/(const double a)const {
		return Body{ this->x / a, this->v / a };
	}
};

void f(std::vector<Body>& body, std::vector<Body>& dot) {
	double k1 = 6, k2 = 6, b1 = 0.1, b2 = 0.2, m1 = 2, m2 = 2;

	int len = 100;

	dot[0].x = body[0].v;
	dot[0].v = (-k1 * (body[0].x - len) - b1 * body[0].v + k2 * (body[1].x - body[0].x - len) + b2 * body[1].v) / m1;
	dot[1].x = body[1].v;
	dot[1].v = (-k2 * (body[1].x - body[0].x - len) - b2 * body[1].v) / m2;
}

void step(std::vector<Body>& st, std::vector<Body>& body, std::vector<Body> k, double h) {
	for (int i = 0; i < st.size(); i++) {
		st[i] = body[i] + k[i] * h / 2;
	}
}

void RungeKutta(std::vector<Body>& body, double h, void (*f)(std::vector<Body>&, std::vector<Body>&)) {
	std::vector<Body> k1, k2, k3, k4, st;
	for (int i = 0; i < body.size(); ++i)
	{
		k1.push_back(Body());
		k2.push_back(Body());
		k3.push_back(Body());
		k4.push_back(Body());
		st.push_back(Body());
	}
	f(body,k1);
	step(st,body,k1,h);
	f(st,k2);

	step(st,body,k2,h);
	f(st,k3);

	step(st,body,k3,h*2);
	f(st,k4);

	for (int i = 0; i < body.size(); ++i) {
		body[i] = body[i] + ((k1[i] + k2[i] * 2 + k3[i] * 2 + k4[i]) * h / 6);
	}
}

int main() {
	sf::RenderWindow win(sf::VideoMode(980, 720), "Window");
	win.setFramerateLimit(60);
	std::ofstream file("time.txt");

	sf::Image spring_image;
	spring_image.loadFromFile("images/spring.png");
	sf::Texture spring_texture;
	spring_texture.loadFromImage(spring_image);
	sf::Sprite spring1, spring2;
	spring1.setTexture(spring_texture);
	spring2.setTexture(spring_texture);
	spring1.setTextureRect(sf::IntRect(0, 0, 502, 188));
	spring2.setTextureRect(sf::IntRect(0, 0, 502, 188));
	spring1.setScale(sf::Vector2f(100.0/502,50.0/188/2));
	spring2.setScale(sf::Vector2f(100.0 / 502, 50.0 / 188/2));
	spring1.setOrigin(502/2, 188/2);
	spring2.setOrigin(502/2, 188/2);

	sf::Vertex lineX[] =
	{
		sf::Vertex(sf::Vector2f(1000,0)),
		sf::Vertex(sf::Vector2f(-1000,0))
	};
	lineX->color= sf::Color(100, 100, 100);
	sf::Vertex lineY[] =
	{
		sf::Vertex(sf::Vector2f(0,-1000)),
		sf::Vertex(sf::Vector2f(0,1000))
	};
	lineY->color = sf::Color(100, 100, 100);

	double h = 0.025;
	double t = 0;

	std::vector<Body> body;
	body.push_back(Body(150, 0));
	body.push_back(Body(350, 0));


	sf::RectangleShape obj1, obj2;
	obj1.setOrigin(25,25);
	obj1.setSize(sf::Vector2f(50, 50));
	obj1.setFillColor(sf::Color::Red);
	obj2.setOrigin(25, 25);
	obj2.setSize(sf::Vector2f(50, 50));
	obj2.setFillColor(sf::Color::Blue);


	sf::View view = win.getDefaultView();
	view.setCenter(300, 0);
	win.setView(view);

	while (win.isOpen()) {
		sf::Event e;
		while (win.pollEvent(e)) {
			if (e.type == sf::Event::Closed)
				win.close();
		}
		win.clear();

		win.draw(lineX, 2, sf::Lines);
		win.draw(lineY, 2, sf::Lines);

		obj1.setPosition(body[0].x, 0);
		obj2.setPosition(body[1].x, 0);
		spring1.setScale((body[0].x - 25) / 502, 50.0 / 188/2);
		spring1.setPosition((body[0].x - 25)/2, 0);
		spring2.setScale((body[1].x - body[0].x -50) / 502, 50.0 / 188/ 2);
		spring2.setPosition((body[1].x + body[0].x) / 2, 0);

		win.draw(obj1);
		win.draw(spring1);
		win.draw(obj2);
		win.draw(spring2);

		RungeKutta(body, h, f);
		t += h;
		if (t - int(t) < h*h)
		{
			file << t << " " << body[0].x << " " << body[0].v << " " << body[1].x << " " << body[1].v << '\n';
		}
		
		win.display();
	}
	file.close();
	return 0;
}
