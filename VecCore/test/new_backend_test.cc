//
// New backend tests
//

#undef NDEBUG
#include <cassert>

#include "backend/scalar.h"
#include "backend/vector.h"

#include <cstdio>

template <class Backend> void test_backend()
{
    typedef typename Backend::Real_t Real_t;
    typedef typename Backend::Real_v Real_v;
    typedef typename Backend::Real_v::Mask Mask_t;

    Mask_t m_true(true), m_false(false);

    assert( m_true.isFull()  == true);
    assert(m_false.isEmpty() == true);

    printf("Real_v::Size       == %d\n\n", Real_v::Size);
    printf("Real_v::Mask::Size == %d\n\n", Real_v::Size);

    Real_v x(0.0);

    for (int i = 0; i < Real_v::Size; i++)
	    x[i] = i;

    Mask_t mask = x > Real_v(Real_v::Size/2.0);

    Real_t *xptr = (Real_t*)(&x);

    for (int i = 0; i < Real_v::Size; i++)
	    printf("x[%d] = %.1f, x[i] > %d == %s\n",
	        i, xptr[i], Real_v::Size/2, mask[i] ? "true" : "false");

    printf("\n");

    for (int i = 0; i < Real_v::Size; i++)
	    x[i] = Real_v::Size - i;

    mask = x > Real_v(Real_v::Size/2.0);

    for (int i = 0; i < Real_v::Size; i++)
	    printf("x[%d] = %.1f, x[i] > %d == %s\n",
	        i, xptr[i], Real_v::Size/2, mask[i] ? "true" : "false");
}

int main(int argc, char *argv[])
{

    printf("Scalar Backend Test (float):\n\n");
    test_backend<VecCore::Backend::Scalar<float>>();

    printf("\n\n");

    printf("Scalar Backend Test (double):\n\n");
    test_backend<VecCore::Backend::Scalar<double>>();

    printf("\n\n");

    printf("Vector Backend Test (float):\n\n");
    test_backend<VecCore::Backend::Vector<float>>();

    printf("\n\n");

    printf("Vector Backend Test (double):\n\n");
    test_backend<VecCore::Backend::Vector<double>>();

    return 0;
}
