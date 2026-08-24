import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';

export default defineConfig({
  plugins: [react()],
  server: {
    port: 5173,
    host: true,
    proxy: {
      '/organisms': 'http://localhost:8000',
      '/report': 'http://localhost:8000'
    }
  }
});
