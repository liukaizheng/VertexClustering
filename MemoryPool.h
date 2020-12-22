#ifndef MEMORYPOOL_H
#define MEMORYPOOL_H

#include <functional>

namespace sigma
{
class PooledAllocator {
    const std::size_t WORDSIZE = 4;
    const std::size_t BLOCKSIZE = 8192;

    std::size_t mRemaining;
    void* mBase;
    void* mLoc;

    void Init()
    {
        mRemaining = 0;
        mBase = nullptr;
        mUsed = 0;
    }
public:
    PooledAllocator() { Init(); }
    ~PooledAllocator() { FreeAll(); }
    void FreeAll()
    {
        while (mBase)
        {
            void *prev = *static_cast<void**>(mBase);
            ::free(mBase);
            mBase = prev;
        }
        Init();
    }

    void* Malloc(const std::size_t req_size)
    {
        const std::size_t size = (req_size + (WORDSIZE - 1)) & ~(WORDSIZE - 1);
        if (size > mRemaining) 
        {
            const std::size_t block_size = 
                (size + sizeof(void *) + (WORDSIZE - 1) > BLOCKSIZE)
                    ? size + sizeof(void *) + (WORDSIZE - 1)
                    : BLOCKSIZE;
            void* m = ::malloc(block_size);
            static_cast<void **>(m)[0] = mBase;
            mBase = m;
            mRemaining = block_size - sizeof(void *);
            mLoc = static_cast<char*>(m) + sizeof(void *);
        }
        void* rloc = mLoc;
        mLoc = static_cast<char*>(mLoc) + size;
        mRemaining -= size;
        mUsed += size;
        return rloc;
    }

    template <typename T>
    T* Allocate(const std::size_t count = 1)
    {
        T* mem = static_cast<T*>(this->Malloc(sizeof(T) * count));
        return mem;
    }

    std::size_t mUsed;
};
}

#endif
