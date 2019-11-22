#include <chrono>

#include "utils.hpp"

int main() {
  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "Testing the Solar Model routines..." << std::endl;

  SolarModel s ("data/SolarModel_AGSS09met_old.dat");
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "Setting up the Solar model took "
            << std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count() << " seconds." << std::endl;

  double test_compton = s.Gamma_P_Compton (3.0, 0.3, 1);
  std::cout << "Test Gamma_Compton value: " << test_compton << std::endl;

  double test_ee = s.Gamma_P_ee (3.0, 0.3);
  std::cout << "Test Gamma_bremsstrahlung value: " << test_ee << std::endl;

  double test_element = s.Gamma_P_element (3.0, 0.3, 1);
  std::cout << "Test element value: " << test_element << std::endl;

  std::cout << "Finished testing!" << std::endl;
  return 0;
}
