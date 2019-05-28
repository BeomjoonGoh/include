// Variadic template
template <typename... Ts> void Print(const Ts& ... args)
{ // loop args using the initialization list, with comma operator.
  double a;
  int dummy[] = { 0,
    ( 
     args.Print(),
     a = args.bar(),
     0)...
  };
  dummy[0]++; // stop complain about unused variable
  std::cout << a << std::endl;
}
