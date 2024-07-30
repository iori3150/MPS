#include <iostream>
#include <yaml-cpp/yaml.h>

int main() {
    // YAMLファイルを読み込む
    YAML::Node config = YAML::LoadFile("config.yaml");
    std::cout <<"config.yaml" << std::endl;

    // キー "name" と "age" を読み込む
    std::string name = config["name"].as<std::string>();
    int age = config["age"].as<int>();

    // 読み込んだ値を表示する
    std::cout << "Name: " << name << std::endl;
    std::cout << "Age: " << age << std::endl;

    return 0;
}
