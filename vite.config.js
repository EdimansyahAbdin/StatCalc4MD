import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

export default defineConfig({
  plugins: [react()],
  // GitHub Pages serves from /StatCalc4MD/ — this must match your repo name exactly
  base: '/StatCalc4MD/',
  build: {
    outDir: 'dist',
    chunkSizeWarningLimit: 900,
  },
  server: {
    port: 3000,
    host: true,
  },
})
