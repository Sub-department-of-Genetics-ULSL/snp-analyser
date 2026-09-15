import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';

export default defineConfig({
  plugins: [react()],
  server: {
    port: 5173,
    host: true,
    proxy: {
      '/api': {
        target: 'http://localhost:8000',
        rewrite: (path) => path.replace(/^\/api/, '')
      },
      '/generated-reports': 'http://localhost:8000',
      '/pdb-files': 'http://localhost:8000'
    }
  }
});
