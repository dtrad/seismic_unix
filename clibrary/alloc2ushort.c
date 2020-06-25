/* allocate a 2-d array */
void **alloc2ushort (size_t n1, size_t n2)
{
	size_t i2;
	void **p;
	size_t size=sizeof(unsigned short);
	if ((p=(void**)malloc(n2*sizeof(void*)))==NULL) 
		return NULL;
	if ((p[0]=(void*)malloc(n2*n1*size))==NULL) {
		free(p);
		return NULL;
	}
	for (i2=0; i2<n2; i2++)
		p[i2] = (char*)p[0]+size*n1*i2;
	return p;
}
void free2ushort(unsigned short **p)
{
	free2((void**)p);
}
